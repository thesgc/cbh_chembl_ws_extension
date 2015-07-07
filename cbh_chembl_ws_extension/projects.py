from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from django.conf import settings
from django.conf.urls import *
from django.http import HttpResponse

from tastypie.resources import ModelResource, Resource
from itertools import chain


from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch, Project, CustomFieldConfig, PinnedCustomField
from cbh_chembl_ws_extension.base import UserResource
from cbh_chembl_ws_extension.authorization import ProjectAuthorization, ProjectListAuthorization
from tastypie.serializers import Serializer
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
import json
import copy
import time
import urllib
from django.core.urlresolvers import reverse

import elasticsearch_client

class ProjectResource(ModelResource):

    class Meta:
        queryset = Project.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get']        
        resource_name = 'cbh_projects'
        authorization = ProjectListAuthorization()
        include_resource_uri = False
        default_format = 'application/json'
        serializer = Serializer()
        filtering = {
            
            "project_key": ALL_WITH_RELATIONS,
        }

    def prepend_urls(self):
        return [
        url(r"^(?P<resource_name>%s)/reindex_elasticsearch/$" % self._meta.resource_name,
                self.wrap_view('reindex_elasticsearch'), name="api_projects_reindex_elasticsearch"),
        ]

    def get_searchform(self, bundle,searchfield_items ):
        '''Note that the form here is expected to have the UOx id as the first item'''
        return {  "cf_form": [
                      {
                          "htmlClass": "col-sm-10",
                          "key": "search_custom_fields__kv_any",
                          #"ngModelOptions": { "updateOn": 'blur' },
                          "disableSuccessState": True,
                          "feedback": False,
                          "options": {
                          "refreshDelay": 0,
                            "async": {
                                "url": reverse("api_get_elasticsearch_autocomplete", 
                                  kwargs={"resource_name": "cbh_compound_batches",
                                  "api_name" : settings.WEBSERVICES_NAME}) ,
                              }
                          }
                      },
                  ],
                  "cf_schema": {
                    "required": [
                                ],
                                "type": "object",
                                "properties": {
                                               "search_custom_fields__kv_any": { 
                                                  "type": "array", 
                                                  "format" : "uiselect",
                                                  "items" :[],
                                                  "placeholder": "Tagged fields",
                                                  "title": "Any of the following custom field values:",
                                                }
                                }
                  },

                  "form": [     
                                
                                {"key": "related_molregno__chembl__chembl_id__in",
                                    "title" : "%s ID" % settings.ID_PREFIX,
                                    "htmlClass": "col-xs-12",
                                    "placeholder" : "Search multiple IDs",
                                    "feedback": False,
                                    "description": "Add your ids and click on them as they appear in the drop-down.",
                                    "options": {
                                      "refreshDelay": 0,
                                      "async": {

                                                    "url": reverse("api_get_elasticsearch_ids", 
                                                      kwargs={"resource_name": "cbh_compound_batches",
                                                      "api_name" : settings.WEBSERVICES_NAME}) ,
                                            }
                                    }},
                                {
                                  "key": "project__project_key__in",
                                  "type": "checkboxes",
                                 # "description": "<info-box freetext='Limit your search results to items tagged with project-specific information'></info-box>",
                                  "placeholder": "Check the console",
                                  "htmlClass": "col-sm-12",
                                  "onChange": "getSearchCustomFields()",
                                  "titleMap": {
                                        p.obj.project_key : p.obj.name for p in  bundle["objects"]
                                  },
                                  
                                  "disableSuccessState": True,
                                },
                                {
                                  "key": "dateStart",
                                  'type': 'datepicker', 
                                  "minDate": "2004-01-01",
                                  "htmlClass": "col-sm-6",
                                  "disableSuccessState": True,
                                  "feedback": False,
                                  'pickadate': {
                                    'selectYears': True, 
                                    'selectMonths': True,
                                  },
                                },
                                {
                                  "key": "dateEnd",
                                  'type': 'datepicker',
                                  "minDate": "2004-01-01",
                                  "htmlClass": "col-sm-6",
                                  "disableSuccessState": True,
                                  "feedback": False,
                                  'pickadate': {
                                    'selectYears': True, 
                                    'selectMonths': True,
                                  },
                                },
                                {
                                    "htmlClass": "col-sm-6",
                                    "disableSuccessState": True,
                                    "feedback": False,
                                    "key": "functional_group",

                                },
                                {
                                  "key": "smiles",
                                  "placeholder": "Search SMILES or SMARTS string",
                                  "append": "today",
                                  "feedback": False,
                                  "htmlClass": "col-sm-6",
                                  "disableSuccessState": True,
                                },                               
                                {
                                  "key": "substruc",
                                  "style": {
                                    "selected": "btn-success",
                                    "unselected": "btn-default"
                                  },
                                  "htmlClass": "col-sm-6",
                                  "type": "radiobuttons",
                                  "disableSuccessState": True,
                                  "feedback": False,
                                  "titleMap": [
                                    {
                                      "value": "with_substructure",
                                      "name": "Substructure",
                                    },
                                    {
                                      "value": "flexmatch",
                                      "name": "Exact Match"
                                    }
                                  ]
                                },
                                {
                                  "key": "multiple_batch_id",
                                  "htmlClass": "col-xs-6",
                                  "disableSuccessState": True,
                                  "feedback": False,

                                },
                                {
                                    "htmlClass": "col-sm-12",
                                    "key": "search_custom_fields__kv_any",
                                    "disableSuccessState": True,
                                    "feedback": False,
                                    "options": {
                                    "refreshDelay": 0,
                                      "async": {
                                          "url": reverse("api_get_elasticsearch_autocomplete", 
                                            kwargs={"resource_name": "cbh_compound_batches",
                                            "api_name" : settings.WEBSERVICES_NAME}) ,
                                        }
                                    }
                                },
                            ],
                            "schema": {
                                "required": [
                                ],
                                "type": "object",
                                "properties": {
                                               "related_molregno__chembl__chembl_id__in":{
                                                  "type": "array", 
                                                  "format" : "uiselect", 
                                                "placeholder": "Look for ids",
                                                "title": "",
                                                "items" :[],
                                                },
                                                "project__project_key__in": {
                                                  "title": "Project",
                                                  "type": "string",
                                                    "enum" :[p.obj.project_key for p in bundle["objects"]],
                                                 "default" : [p.obj.project_key for p in bundle["objects"]]

                                                },
                                                "multiple_batch_id": {
                                                  "title": "Multiple Batch ID",
                                                  "type": "string",
                                                },
                                                "functional_group": {
                                                    "title" : "Functional Group",
                                                    "type" : "string",
                                                     "format" : "uiselect", 
                                                      "placeholder": "Search chemical groups",
                                                      "options": {
                                                          "searchDescriptions": False
                                                      },
                                                      "default" : "",
                                                      "copyValueTo" : "smiles",
                                                    "items" : [{
                                                                    "label":"None",
                                                                    "value":"",
                                                                   # "description":""
                                                                  },] + sorted([
                                                                  
                                                                  {
                                                                    "label":"Alkyl Carbon",
                                                                    "value":"[CX4]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Allenic Carbon",
                                                                    "value":"[$([CX2](=C)=C)]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Vinylic Carbon",
                                                                    "value":"[$([CX3]=[CX3])]",
                                                                   # "description":"Ethenyl carbon"
                                                                  },
                                                                  {
                                                                    "label":"Acetylenic Carbon",
                                                                    "value":"[$([CX2]#C)]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Arene",
                                                                    "value":"c",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Carbonyl group. Low specificity",
                                                                    "value":"[CX3]=[OX1]",
                                                                   # "description":"Hits carboxylic acid, ester, ketone, aldehyde, carbonic acid/ester,anhydride, carbamic acid/ester, acyl halide, amide."
                                                                  },
                                                                  {
                                                                    "label":"Carbonyl group",
                                                                    "value":"[$([CX3]=[OX1]),$([CX3+]-[OX1-])]",
                                                                   # "description":"Hits either resonance structure"
                                                                  },
                                                                  {
                                                                    "label":"Carbonyl with Carbon",
                                                                    "value":"[CX3](=[OX1])C",
                                                                   # "description":"Hits aldehyde, ketone, carboxylic acid (except formic), anhydride (except formic), acyl halides (acid halides). Won't hit carbamic acid/ester, carbonic acid/ester."
                                                                  },
                                                                  {
                                                                    "label":"Carbonyl with Nitrogen.",
                                                                    "value":"[OX1]=CN",
                                                                   # "description":"Hits amide, carbamic acid/ester, poly peptide"
                                                                  },
                                                                  {
                                                                    "label":"Carbonyl with Oxygen.",
                                                                    "value":"[CX3](=[OX1])O",
                                                                   # "description":"Hits ester, carboxylic acid, carbonic acid or ester, carbamic acid or ester, anhydride Won't hit aldehyde or ketone."
                                                                  },
                                                                  {
                                                                    "label":"Acyl Halide",
                                                                    "value":"[CX3](=[OX1])[F,Cl,Br,I]",
                                                                   # "description":"acid halide, -oyl halide"
                                                                  },
                                                                  {
                                                                    "label":"Aldehyde",
                                                                    "value":"[CX3H1](=O)[#6]",
                                                                   # "description":"-al"
                                                                  },
                                                                  {
                                                                    "label":"Anhydride",
                                                                    "value":"[CX3](=[OX1])[OX2][CX3](=[OX1])",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Amide",
                                                                    "value":"[NX3][CX3](=[OX1])[#6]",
                                                                   # "description":"-amide"
                                                                  },
                                                                  {
                                                                    "label":"Amidinium",
                                                                    "value":"[NX3][CX3]=[NX3+]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Carbamate.",
                                                                    "value":"[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]",
                                                                   # "description":"Hits carbamic esters, acids, and zwitterions"
                                                                  },
                                                                  {
                                                                    "label":"Carbamic ester",
                                                                    "value":"[NX3][CX3](=[OX1])[OX2H0]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Carbamic acid.",
                                                                    "value":"[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]",
                                                                   # "description":"Hits carbamic acids and zwitterions."
                                                                  },
                                                                  {
                                                                    "label":"Carboxylate Ion.",
                                                                    "value":"[CX3](=O)[O-]",
                                                                   # "description":"Hits conjugate bases of carboxylic, carbamic, and carbonic acids."
                                                                  },
                                                                  {
                                                                    "label":"Carbonic Acid or Carbonic Ester",
                                                                    "value":"[CX3](=[OX1])(O)O",
                                                                   # "description":"Carbonic Acid, Carbonic Ester, or combination"
                                                                  },
                                                                  {
                                                                    "label":"Carbonic Acid or Carbonic Acid-Ester",
                                                                    "value":"[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]",
                                                                   # "description":"Hits acid and conjugate base. Won't hit carbonic acid diester"
                                                                  },
                                                                  {
                                                                    "label":"Carbonic Ester (carbonic acid diester)",
                                                                    "value":"C[OX2][CX3](=[OX1])[OX2]C",
                                                                   # "description":"Won't hit carbonic acid or combination carbonic acid/ester"
                                                                  },
                                                                  {
                                                                    "label":"Carboxylic acid",
                                                                    "value":"[CX3](=O)[OX2H1]",
                                                                   # "description":"-oic acid, COOH"
                                                                  },
                                                                  {
                                                                    "label":"Carboxylic acid or conjugate base.",
                                                                    "value":"[CX3](=O)[OX1H0-,OX2H1]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Cyanamide",
                                                                    "value":"[NX3][CX2]#[NX1]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Ester Also hits anhydrides",
                                                                    "value":"[#6][CX3](=O)[OX2H0][#6]",
                                                                   # "description":"won't hit formic anhydride."
                                                                  },
                                                                  {
                                                                    "label":"Ketone",
                                                                    "value":"[#6][CX3](=O)[#6]",
                                                                   # "description":"-one"
                                                                  },
                                                                  {
                                                                    "label":"Ether",
                                                                    "value":"[OD2]([#6])[#6]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Hydrogen Atom",
                                                                    "value":"[H]",
                                                                   # "description":"Hits SMILES that are hydrogen atoms: [H+] [2H] [H][H]"
                                                                  },
                                                                  {
                                                                    "label":"Not a Hydrogen Atom",
                                                                    "value":"[!#1]",
                                                                   # "description":"Hits SMILES that are not hydrogen atoms."
                                                                  },
                                                                  {
                                                                    "label":"Proton",
                                                                    "value":"[H+]",
                                                                   # "description":"Hits positively charged hydrogen atoms: [H+]"
                                                                  },
                                                                  {
                                                                    "label":"Mono-Hydrogenated Cation",
                                                                    "value":"[+H]",
                                                                   # "description":"Hits atoms that have a positive charge and exactly one attached hydrogen: F[C+](F)[H]"
                                                                  },
                                                                  {
                                                                    "label":"Not Mono-Hydrogenated",
                                                                    "value":"[!H] or [!H1]",
                                                                   # "description":"Hits atoms that don't have exactly one attached hydrogen.</p></dl><a NAME=\"N\"></a> see carbonyl</p>mine (-amino) </h3><dl>"
                                                                  },
                                                                  {
                                                                    "label":"Primary or secondary amine, not amide.",
                                                                    "value":"[NX3;H2,H1;!$(NC=O)]",
                                                                   # "description":"Not ammonium ion (N must be 3-connected), not ammonia (H count can't be 3). Primary or secondary is specified by N's H-count (H2 & H1 respectively). Also note that \"&\" (and) is the dafault opperator and is higher precedence that \",\" (or), which is higher precedence than \";\" (and). Will hit cyanamides and thioamides"
                                                                  },
                                                                  {
                                                                    "label":"Enamine",
                                                                    "value":"[NX3][CX3]=[CX3]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Primary amine, not amide.",
                                                                    "value":"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6] Not amide (C not double bonded to a hetero-atom), not ammonium ion (N must be 3-connected), not ammonia (N's H-count can't be 3), not cyanamide (C not triple bonded to a hetero-atom)",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Two primary or secondary amines",
                                                                    "value":"[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]",
                                                                   # "description":"Here we use the disconnection symbol (\".\") to match two separate unbonded identical patterns."
                                                                  },
                                                                  {
                                                                    "label":"Enamine or Aniline Nitrogen",
                                                                    "value":"[NX3][$(C=C),$(cc)]",
                                                                   # "description":""
                                                                  },
                                                                  # {
                                                                  #   "label":"Generic amino acid: low specificity.",
                                                                  #   "value":"[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]",
                                                                  #  # "description":"For use w/ non-standard a.a. search. hits pro but not gly. Hits acids and conjugate bases. Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Dipeptide group. generic amino acid: low specificity.",
                                                                  #   "value":"[NX3H2,NH3X4+][CX4H]([*])[CX3](=[OX1])[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-]",
                                                                  #  # "description":"Won't hit pro or gly. Hits acids and conjugate bases."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Amino Acid",
                                                                  #   "value":"[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]",
                                                                  #  # "description":"Replace * w/ a specific a.a. side chain from the 18_standard_side_chains list to hit a specific standard a.a. Won't work with Proline or Glycine, they have their own SMARTS (see side chain list). Hits acids and conjugate bases. Hits single a.a.s and specific residues w/i n polypeptides (internal, or terminal). {e.g. usage: Alanine side chain is [CH3X4] . Search is [$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([ CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]}"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Alanine side chain",
                                                                  #   "value":"[CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Arginine side chain.",
                                                                  #   "value":"[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]",
                                                                  #  # "description":"Hits acid and conjugate base."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Aspargine side chain.",
                                                                  #   "value":"[CH2X4][CX3](=[OX1])[NX3H2]",
                                                                  #  # "description":"Also hits Gln side chain when used alone."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Aspartate (or Aspartic acid) side chain.",
                                                                  #   "value":"[CH2X4][CX3](=[OX1])[OH0-,OH]",
                                                                  #  # "description":"Hits acid and conjugate base. Also hits Glu side chain when used alone."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Cysteine side chain.",
                                                                  #   "value":"[CH2X4][SX2H,SX1H0-]",
                                                                  #  # "description":"Hits acid and conjugate base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Glutamate (or Glutamic acid) side chain.",
                                                                  #   "value":"[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]",
                                                                  #  # "description":"Hits acid and conjugate base."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Glycine",
                                                                  #   "value":"[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Histidine side chain.",
                                                                  #   "value":"[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:<br>[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1",
                                                                  #  # "description":"Hits acid & conjugate base for either Nitrogen. Note that the Ns can be either ([(Cationic 3-connected with one H) or (Neutral 2-connected without any Hs)] where there is a second-neighbor who is [3-connected with one H]) or (3-connected with one H)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Isoleucine side chain",
                                                                  #   "value":"[CHX4]([CH3X4])[CH2X4][CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Leucine side chain",
                                                                  #   "value":"[CH2X4][CHX4]([CH3X4])[CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Lysine side chain.",
                                                                  #   "value":"[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]",
                                                                  #  # "description":"Acid and conjugate base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Methionine side chain",
                                                                  #   "value":"[CH2X4][CH2X4][SX2][CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Phenylalanine side chain",
                                                                  #   "value":"[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Proline",
                                                                  #   "value":"[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Serine side chain",
                                                                  #   "value":"[CH2X4][OX2H]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Thioamide",
                                                                  #   "value":"[NX3][CX3]=[SX1]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Threonine side chain",
                                                                  #   "value":"[CHX4]([CH3X4])[OX2H]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Tryptophan side chain",
                                                                  #   "value":"[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Tyrosine side chain.",
                                                                  #   "value":"[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1",
                                                                  #  # "description":"Acid and conjugate base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Valine side chain",
                                                                  #   "value":"[CHX4]([CH3X4])[CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Alanine side chain",
                                                                  #   "value":"[CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Arginine side chain.",
                                                                  #   "value":"[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]",
                                                                  #  # "description":"Hits acid and conjugate base."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Aspargine side chain.",
                                                                  #   "value":"[CH2X4][CX3](=[OX1])[NX3H2]",
                                                                  #  # "description":"Also hits Gln side chain when used alone."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Aspartate (or Aspartic acid) side chain.",
                                                                  #   "value":"[CH2X4][CX3](=[OX1])[OH0-,OH]",
                                                                  #  # "description":"Hits acid and conjugate base. Also hits Glu side chain when used alone."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Cysteine side chain.",
                                                                  #   "value":"[CH2X4][SX2H,SX1H0-]",
                                                                  #  # "description":"Hits acid and conjugate base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Glutamate (or Glutamic acid) side chain.",
                                                                  #   "value":"[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]",
                                                                  #  # "description":"Hits acid and conjugate base."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Glycine",
                                                                  #   "value":"N[CX4H2][CX3](=[OX1])[O,N]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Histidine side chain.",
                                                                  #   "value":"[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:<br>[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1",
                                                                  #  # "description":"Hits acid & conjugate base for either Nitrogen. Note that the Ns can be either ([(Cationic 3-connected with one H) or (Neutral 2-connected without any Hs)] where there is a second-neighbor who is [3-connected"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Isoleucine side chain",
                                                                  #   "value":"[CHX4]([CH3X4])[CH2X4][CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Leucine side chain",
                                                                  #   "value":"[CH2X4][CHX4]([CH3X4])[CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Lysine side chain.",
                                                                  #   "value":"[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]",
                                                                  #  # "description":"Acid and conjugate base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Methionine side chain",
                                                                  #   "value":"[CH2X4][CH2X4][SX2][CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Phenylalanine side chain",
                                                                  #   "value":"[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Proline",
                                                                  #   "value":"N1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[O,N]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Serine side chain",
                                                                  #   "value":"[CH2X4][OX2H]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Threonine side chain",
                                                                  #   "value":"[CHX4]([CH3X4])[OX2H]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Tryptophan side chain",
                                                                  #   "value":"[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Tyrosine side chain.",
                                                                  #   "value":"[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1",
                                                                  #  # "description":"Acid and conjugate base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Valine side chain",
                                                                  #   "value":"[CHX4]([CH3X4])[CH3X4]",
                                                                  #  # "description":""
                                                                  # },
                                                                  {
                                                                    "label":"Azide group.",
                                                                    "value":"[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]",
                                                                   # "description":"Hits any atom with an attached azide."
                                                                  },
                                                                  {
                                                                    "label":"Azide ion.",
                                                                    "value":"[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]",
                                                                   # "description":"Hits N in azide ion"
                                                                  },
                                                                  {
                                                                    "label":"Nitrogen.",
                                                                    "value":"[#7]",
                                                                   # "description":"Nitrogen in N-containing compound. aromatic or aliphatic. Most general interpretation of \"azo\""
                                                                  },
                                                                  {
                                                                    "label":"Azo Nitrogen. Low specificity.",
                                                                    "value":"[NX2]=N",
                                                                   # "description":"Hits diazene, azoxy and some diazo structures"
                                                                  },
                                                                  {
                                                                    "label":"Azo Nitrogen.diazene",
                                                                    "value":"[NX2]=[NX2]",
                                                                   # "description":"(diaza alkene)"
                                                                  },
                                                                  {
                                                                    "label":"Azoxy Nitrogen.",
                                                                    "value":"[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Diazo Nitrogen",
                                                                    "value":"[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Azole.",
                                                                    "value":"[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]",
                                                                   # "description":"5 member aromatic heterocycle w/ 2double bonds. contains N & another non C (N,O,S) subclasses are furo-, thio-, pyrro- (replace CH o' furfuran, thiophene, pyrrol w/ N)"
                                                                  },
                                                                  {
                                                                    "label":"Hydrazine H2NNH2",
                                                                    "value":"[NX3][NX3]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Hydrazone C=NNH2",
                                                                    "value":"[NX3][NX2]=[*]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Substituted imine",
                                                                    "value":"[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]",
                                                                   # "description":"Schiff base"
                                                                  },
                                                                  {
                                                                    "label":"Substituted or un-substituted imine",
                                                                    "value":"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Iminium",
                                                                    "value":"[NX3+]=[CX3]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Unsubstituted dicarboximide",
                                                                    "value":"[CX3](=[OX1])[NX3H][CX3](=[OX1])",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Substituted dicarboximide",
                                                                    "value":"[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Dicarboxdiimide",
                                                                    "value":"[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Nitrate group",
                                                                    "value":"[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]",
                                                                   # "description":"Also hits nitrate anion"
                                                                  },
                                                                  {
                                                                    "label":"Nitrate Anion",
                                                                    "value":"[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Nitrile",
                                                                    "value":"[NX1]#[CX2]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Isonitrile",
                                                                    "value":"[CX1-]#[NX2+]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Nitro group.",
                                                                    "value":"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8] Hits both forms.",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Two Nitro groups",
                                                                    "value":"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8].[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Nitroso-group",
                                                                    "value":"[NX2]=[OX1]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"N-Oxide",
                                                                    "value":"[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
                                                                   # "description":"Hits both forms. Won't hit azoxy, nitro, nitroso,or nitrate."
                                                                  },
                                                                  {
                                                                    "label":"Hydroxyl",
                                                                    "value":"[OX2H]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Hydroxyl in Alcohol",
                                                                    "value":"[#6][OX2H]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Hydroxyl in Carboxylic Acid",
                                                                    "value":"[OX2H][CX3]=[OX1]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Hydroxyl in H-O-P-",
                                                                    "value":"[OX2H]P",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Enol",
                                                                    "value":"[OX2H][#6X3]=[#6]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Phenol",
                                                                    "value":"[OX2H][cX3]:[c]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Enol or Phenol",
                                                                    "value":"[OX2H][$(C=C),$(cc)]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Hydroxyl_acidic",
                                                                    "value":"[$([OH]-*=[!#6])]",
                                                                   # "description":"An acidic hydroxyl is a hydroxyl bonded to an atom which is multiply bonded to a hetero atom, this includes carboxylic, sulphur, phosphorous, halogen and nitrogen oxyacids."
                                                                  },
                                                                  {
                                                                    "label":"Peroxide groups.",
                                                                    "value":"[OX2,OX1-][OX2,OX1-]",
                                                                   # "description":"Also hits anions."
                                                                  },
                                                                  {
                                                                    "label":"Phosphoric_acid groups.",
                                                                    "value":"[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
                                                                   # "description":"Hits both depiction forms. Hits orthophosphoric acid and polyphosphoric acid anhydrides. Doesn't hit monophosphoric acid anhydride esters (including acidic mono- & di- esters) but will hit some polyphosphoric acid anhydride esters (mono- esters on pyrophosphoric acid and longer, di- esters on linear triphosphoric acid and longer)."
                                                                  },
                                                                  {
                                                                    "label":"Phosphoric_ester groups.",
                                                                    "value":"[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
                                                                   # "description":"Hits both depiction forms. Doesn't hit non-ester phosphoric_acid groups."
                                                                  },
                                                                  {
                                                                    "label":"Carbo-Thiocarboxylate",
                                                                    "value":"[S-][CX3](=S)[#6]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Carbo-Thioester",
                                                                    "value":"S([#6])[CX3](=O)[#6]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Thio analog of carbonyl",
                                                                    "value":"[#6X3](=[SX1])([!N])[!N]",
                                                                   # "description":"Where S replaces O. Not a thioamide."
                                                                  },
                                                                  {
                                                                    "label":"Thiol, Sulfide or Disulfide Sulfur",
                                                                    "value":"[SX2]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Thiol",
                                                                    "value":"[#16X2H]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Sulfur with at-least one hydrogen.",
                                                                    "value":"[#16!H0]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Thioamide",
                                                                    "value":"[NX3][CX3]=[SX1]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Sulfide",
                                                                    "value":"[#16X2H0]",
                                                                   # "description":"-alkylthio Won't hit thiols. Hits disulfides."
                                                                  },
                                                                  {
                                                                    "label":"Mono-sulfide",
                                                                    "value":"[#16X2H0][!#16]",
                                                                   # "description":"alkylthio- or alkoxy- Won't hit thiols. Won't hit disulfides."
                                                                  },
                                                                  {
                                                                    "label":"Di-sulfide",
                                                                    "value":"[#16X2H0][#16X2H0]",
                                                                   # "description":"Won't hit thiols. Won't hit mono-sulfides."
                                                                  },
                                                                  {
                                                                    "label":"Two Sulfides",
                                                                    "value":"[#16X2H0][!#16].[#16X2H0][!#16]",
                                                                   # "description":"Won't hit thiols. Won't hit mono-sulfides. Won't hit disulfides."
                                                                  },
                                                                  {
                                                                    "label":"Sulfinate",
                                                                    "value":"[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]",
                                                                   # "description":"Won't hit Sulfinic Acid. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfinic Acid",
                                                                    "value":"[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]",
                                                                   # "description":"Won't hit substituted Sulfinates. Hits Both Depiction Forms. Hits acid and conjugate base (sulfinate)."
                                                                  },
                                                                  {
                                                                    "label":"Sulfone. Low specificity.",
                                                                    "value":"[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]",
                                                                   # "description":"Hits all sulfones, including heteroatom-substituted sulfones: sulfonic acid, sulfonate, sulfuric acid mono- & di- esters, sulfamic acid, sulfamate, sulfonamide... Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfone. High specificity.",
                                                                    "value":"[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]",
                                                                   # "description":"Only hits carbo- sulfones (Won't hit herteroatom-substituted molecules). Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfonic acid. High specificity.",
                                                                    "value":"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
                                                                   # "description":"Only hits carbo- sulfonic acids (Won't hit herteroatom-substituted molecules). Hits acid and conjugate base. Hits Both Depiction Forms. Hits Arene sulfonic acids."
                                                                  },
                                                                  {
                                                                    "label":"Sulfonate",
                                                                    "value":"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]",
                                                                   # "description":"(sulfonic ester) Only hits carbon-substituted sulfur (Oxygen may be herteroatom-substituted). Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfonamide.",
                                                                    "value":"[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]",
                                                                   # "description":"Only hits carbo- sulfonamide. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Carbo-azosulfone",
                                                                    "value":"[SX4](C)(C)(=O)=N",
                                                                   # "description":"Partial N-Analog of Sulfone"
                                                                  },
                                                                  {
                                                                    "label":"Sulfonamide",
                                                                    "value":"[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]",
                                                                   # "description":"(sulf drugs) Won't hit sulfamic acid or sulfamate. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfoxide Low specificity.",
                                                                    "value":"[$([#16X3]=[OX1]),$([#16X3+][OX1-])]",
                                                                   # "description":"( sulfinyl, thionyl ) Analog of carbonyl where S replaces C. Hits all sulfoxides, including heteroatom-substituted sulfoxides, dialkylsulfoxides carbo-sulfoxides, sulfinate, sulfinic acids... Hits Both Depiction Forms. Won't hit sulfones."
                                                                  },
                                                                  {
                                                                    "label":"Sulfoxide High specificity",
                                                                    "value":"[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]",
                                                                   # "description":"(sulfinyl , thionyl) Analog of carbonyl where S replaces C. Only hits carbo-sulfoxides (Won't hit herteroatom-substituted molecules). Hits Both Depiction Forms. Won't hit sulfones."
                                                                  },
                                                                  {
                                                                    "label":"Sulfate",
                                                                    "value":"[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
                                                                   # "description":"(sulfuric acid monoester) Only hits when oxygen is carbon-substituted. Hits acid and conjugate base. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfuric acid ester (sulfate ester) Low specificity.",
                                                                    "value":"[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]",
                                                                   # "description":"Hits sulfuric acid, sulfuric acid monoesters (sulfuric acids) and diesters (sulfates). Hits acid and conjugate base. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfuric Acid Diester.",
                                                                    "value":"[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]",
                                                                   # "description":"Only hits when oxygen is carbon-substituted. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfamate.",
                                                                    "value":"[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]",
                                                                   # "description":"Only hits when oxygen is carbon-substituted. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfamic Acid.",
                                                                    "value":"[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]",
                                                                   # "description":"Hits acid and conjugate base. Hits Both Depiction Forms."
                                                                  },
                                                                  {
                                                                    "label":"Sulfenic acid.",
                                                                    "value":"[#16X2][OX2H,OX1H0-]",
                                                                   # "description":"Hits acid and conjugate base."
                                                                  },
                                                                  {
                                                                    "label":"Sulfenate.",
                                                                    "value":"[#16X2][OX2H0]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Any carbon attached to any halogen",
                                                                    "value":"[#6][F,Cl,Br,I]",
                                                                   # "description":""
                                                                  },
                                                                  {
                                                                    "label":"Halogen",
                                                                    "value":"[F,Cl,Br,I]",
                                                                   # "description":""
                                                                  },
                                                                  # {
                                                                  #   "label":"Three_halides groups",
                                                                  #   "value":"[F,Cl,Br,I].[F,Cl,Br,I].[F,Cl,Br,I]",
                                                                  #  # "description":"Hits SMILES that have three halides."
                                                                  # },
                                                                  
                                                                  # {
                                                                  #   "label":"sp2 cationic carbon",
                                                                  #   "value":"[$([cX2+](:*):*)]",
                                                                  #  # "description":"Aromatic cationic sp2 carbon with a free electron in a non-bonding sp2 hybrid orbital"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Aromatic sp2 carbon.",
                                                                  #   "value":"[$([cX3](:*):*),$([cX2+](:*):*)]",
                                                                  #  # "description":"The first recursive SMARTS matches carbons that are three-connected, the second case matches two-connected carbons (i.e cations with a free electron in a non-bonding sp2 hybrid orbital)"
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp2 carbon.",
                                                                  #   "value":"[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*)]",
                                                                  #  # "description":"The first recursive SMARTS matches carbons that are three-connected and aromatic. The second case matches two-connected aromatic ca rbons (i.e cations with a free electron in a non-bonding sp2 hybrid orbital). The third case matches three-connected non-aromatic carbons ( alkenes). The fourth case matches non-aromatic cationic alkene carbons."
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp2 nitrogen.",
                                                                  #   "value":"[$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)]",
                                                                  #  # "description":"Can be aromatic 3-connected with 2 aromatic bonds (eg pyrrole,Pyridine-N-oxide), aromatic 2-connected with 2 aromatic bonds (and a free pair of electrons in a nonbonding orbital, e.g.Pyridine), either aromatic or non-aromatic 2-connected with a double bond (and a free pair of electrons in a nonbonding orbital, e.g. C=N ), non aromatic 3-connected with 2 double bonds (e.g. a nitro group; this form does not exist in reality, SMILES can represent the charge-separated resonance structures as a single uncharged structure), either aromatic or non-aromatic 3-connected cation w/ 1 single bond and 1 double bond (e.g. a nitro group, here the individual charge-separated resonance structures are specified), either aromatic or non-aromatic 3-connected hydrogenated cation with a double bond (as the previous case but R is hydrogen), rspectively."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Explicit Hydrogen on sp2-Nitrogen",
                                                                  #   "value":"[$([#1X1][$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)])]",
                                                                  #  # "description":"(H must be an isotope or ion)"
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp3 nitrogen",
                                                                  #   "value":"[$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)]",
                                                                  #  # "description":"One atom that is (a 4-connected N cation or a 3-connected N) and is not double bonded and is not aromatically bonded."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Explicit Hydrogen on an sp3 N.",
                                                                  #   "value":"[$([#1X1][$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)])]",
                                                                  #  # "description":"One atom that is a 1-connected H that is bonded to an sp3 N. (H must be an isotope or ion)"
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp2 N in N-Oxide",
                                                                  #   "value":"[$([$([NX3]=O),$([NX3+][O-])])]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp3 N in N-Oxide Exclusive:",
                                                                  #   "value":"[$([$([NX4]=O),$([NX4+][O-])])]",
                                                                  #  # "description":"Only hits if O is explicitly present. Won't hit if * is in SMILES in place of O."
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp3 N in N-Oxide Inclusive:",
                                                                  #   "value":"[$([$([NX4]=O),$([NX4+][O-,#0])])]",
                                                                  #  # "description":"Hits if O could be present. Hits if * if used in place of O in smiles."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Quaternary Nitrogen",
                                                                  #   "value":"[$([NX4+]),$([NX4]=*)]",
                                                                  #  # "description":"Hits non-aromatic Ns."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Tricoordinate S double bonded to N.",
                                                                  #   "value":"[$([SX3]=N)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"S double-bonded to Carbon",
                                                                  #   "value":"[$([SX1]=[#6])]",
                                                                  #  # "description":"Hits terminal (1-connected S)"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Triply bonded N",
                                                                  #   "value":"[$([NX1]#*)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Divalent Oxygen",
                                                                  #   "value":"[$([OX2])]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Unbranched_alkane groups.",
                                                                  #   "value":"[R0;D2][R0;D2][R0;D2][R0;D2]",
                                                                  #  # "description":"Only hits alkanes (single-bond chains). Only hits chains of at-least 4 members. All non-(implicit-hydrogen) atoms count as branches (e.g. halide substituted chains count as branched)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Unbranched_chain groups.",
                                                                  #   "value":"[R0;D2]~[R0;D2]~[R0;D2]~[R0;D2]",
                                                                  #  # "description":"Hits any bond (single, double, triple). Only hits chains of at-least 4 members. All non-(implicit-hydrogen) atoms count as branches (e.g. halide substituted chains count as branched)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Long_chain groups.",
                                                                  #   "value":"[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]",
                                                                  #  # "description":"Aliphatic chains at-least 8 members long."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Atom_fragment",
                                                                  #   "value":"[!$([#6+0]);!$(C(F)(F)F);!$(c(:[!c]):[!c])!$([#6]=,#[!#6])]",
                                                                  #  # "description":"(CLOGP definition) A fragment atom is a not an isolating carbon"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Carbon_isolating",
                                                                  #   "value":"[$([#6+0]);!$(C(F)(F)F);!$(c(:[!c]):[!c])!$([#6]=,#[!#6])]",
                                                                  #  # "description":"This definition is based on that in CLOGP, so it is a charge-neutral carbon, which is not a CF3 or an aromatic C between two aromati c hetero atoms eg in tetrazole, it is not multiply bonded to a hetero atom."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Terminal S bonded to P",
                                                                  #   "value":"[$([SX1]~P)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Nitrogen on -N-C=N-",
                                                                  #   "value":"[$([NX3]C=N)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Nitrogen on -N-N=C-",
                                                                  #   "value":"[$([NX3]N=C)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Nitrogen on -N-N=N-",
                                                                  #   "value":"[$([NX3]N=N)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Oxygen in -O-C=N-",
                                                                  #   "value":"[$([OX2]C=N)]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Rotatable bond",
                                                                  #   "value":"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",
                                                                  #  # "description":"An atom which is not triply bonded and not one-connected i.e.terminal connected by a single non-ring bond to and equivalent atom. Note that logical operators can be applied to bonds (\"-&!@\"). Here, the overall SMARTS consists of two atoms and one bond. The bond is \"site and not ring\". *#* any atom triple bonded to any atom. By enclosing this SMARTS in parentheses and preceding with $, this enables us to use $(*#*) to write a recursive SMARTS using that string as an atom primitive. The purpose is to avoid bonds such as c1ccccc1-C#C which wo be considered rotatable without this specification."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Bicyclic",
                                                                  #   "value":"[$([*R2]([*R])([*R])([*R]))].[$([*R2]([*R])([*R])([*R]))]",
                                                                  #  # "description":"Bicyclic compounds have 2 bridgehead atoms with 3 arms connecting the bridgehead atoms."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Ortho",
                                                                  #   "value":"*-!:aa-!:*",
                                                                  #  # "description":"Ortho-substituted ring"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Meta",
                                                                  #   "value":"*-!:aaa-!:*",
                                                                  #  # "description":"Meta-substituted ring"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Para",
                                                                  #   "value":"*-!:aaaa-!:*",
                                                                  #  # "description":"Para-substituted ring"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Acylic-bonds",
                                                                  #   "value":"*!@*",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Single bond and not in a ring",
                                                                  #   "value":"*-!@*",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Non-ring atom",
                                                                  #   "value":"[R0] or [!R]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Macrocycle groups.",
                                                                  #   "value":"[r;!r3;!r4;!r5;!r6;!r7]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"S in aromatic 5-ring with lone pair",
                                                                  #   "value":"[sX2r5]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Aromatic 5-Ring O with Lone Pair",
                                                                  #   "value":"[oX2r5]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"N in 5-sided aromatic ring",
                                                                  #   "value":"[nX2r5]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Spiro-ring center",
                                                                  #   "value":"[X4;R2;r4,r5,r6](@[r4,r5,r6])(@[r4,r5,r6])(@[r4,r5,r6])@[r4,r5,r6]rings size 4-6",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"N in 5-ring arom",
                                                                  #   "value":"[$([nX2r5]:[a-]),$([nX2r5]:[a]:[a-])] anion",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"CIS or TRANS double bond in a ring",
                                                                  #   "value":"*/,\\[R]=;@[R]/,\\*",
                                                                  #  # "description":"An isomeric SMARTS consisting of four atoms and three bonds."
                                                                  # },
                                                                  # {
                                                                  #   "label":"CIS or TRANS double or aromatic bond in a ring",
                                                                  #   "value":"*/,\\[R]=,:;@[R]/,\\*",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Unfused benzene ring",
                                                                  #   "value":"[cR1]1[cR1][cR1][cR1][cR1][cR1]1",
                                                                  #  # "description":"To find a benzene ring which is not fused, we write a SMARTS of 6 aromatic carbons in a ring where each atom is only in one ring:"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Multiple non-fused benzene rings",
                                                                  #   "value":"[cR1]1[cR1][cR1][cR1][cR1][cR1]1.[cR1]1[cR1][cR1][cR1][cR1][cR1]1",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Fused benzene rings",
                                                                  #   "value":"c12ccccc1cccc2",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Generic amino acid: low specificity.",
                                                                  #   "value":"[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]",
                                                                  #  # "description":"For use w/ non-standard a.a. search. hits pro but not gly. Hits acids and conjugate bases. Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"A.A. Template for 20 standard a.a.s",
                                                                  #   "value":"[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),<br>$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N])]",
                                                                  #  # "description":"Pro, Gly, Other. Replace * w/ the entire 18_standard_side_chains list to get \"any standard a.a.\" Hits acids and conjugate bases. Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Proline",
                                                                  #   "value":"[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Glycine",
                                                                  #   "value":"[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Other a.a.",
                                                                  #   "value":"[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]",
                                                                  #  # "description":"Replace * w/ a specific a.a. side chain from the 18_standard_side_chains list to hit a specific standard a.a. Won't work with Proline or Glycine, they have their own SMARTS (see side chain list). Hits acids and conjugate bases. Hits single a.a.s and specific residues w/i polypeptides (internal, or terminal). Example usage: Alanine side chain is [CH3X4] Alanine Search is [$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]"
                                                                  # },
                                                                  # {
                                                                  #   "label":"18_standard_aa_side_chains.",
                                                                  #   "value":"([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])",
                                                                  #  # "description":"Can be any of the standard 18 (Pro & Gly are treated separately) Hits acids and conjugate bases."
                                                                  # },
                                                                  # {
                                                                  #   "label":"N in Any_standard_amino_acid.",
                                                                  #   "value":"[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]",
                                                                  #  # "description":"Format is A.A.Template for 20 standard a.a.s. where * is replaced by the entire 18_standard_side_chains list (or'd together). A generic amino acid with any of the 18 side chains or, proline or glycine. Hits \"standard\" amino acids that have terminally appended groups (i.e. \"standard\" refers to the side chains). (Pro, Gly, or 18 normal a.a.s.) Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Non-standard amino acid.",
                                                                  #   "value":"[$([NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]",
                                                                  #  # "description":"Generic amino acid but not a \"standard\" amino acid (\"standard\" refers to the 20 normal side chains). Won't hit amino acids that are non-standard due solely to the fact that groups are terminally-appended to the polypeptide chain (N or C term). format is [$(generic a.a.);!$(not a standard one)] Hits single a.a.s and specific residues w/in polypeptides (internal, or terminal)."
                                                                  # },
                                                                  {
                                                                    "label":"Sulfide",
                                                                    "value":"[#16X2H0]",
                                                                   # "description":"(-alkylthio) Won't hit thiols. Hits disulfides too."
                                                                  },
                                                                  {
                                                                    "label":"Mono-sulfide",
                                                                    "value":"[#16X2H0][!#16]",
                                                                   # "description":"(alkylthio- or alkoxy-) R-S-R Won't hit thiols. Won't hit disulfides."
                                                                  },
                                                                  {
                                                                    "label":"Di-sulfide",
                                                                    "value":"[#16X2H0][#16X2H0]",
                                                                   # "description":"Won't hit thiols. Won't hit mono-sulfides."
                                                                  },
                                                                  # {
                                                                  #   "label":"Two sulfides",
                                                                  #   "value":"[#16X2H0][!#16].[#16X2H0][!#16]",
                                                                  #  # "description":"Won't hit thiols. Won't hit mono-sulfides. Won't hit disulfides."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Acid/conj-base",
                                                                  #   "value":"[OX2H,OX1H0-]",
                                                                  #  # "description":"Hits acid and conjugate base. acid/base"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Non-acid Oxygen",
                                                                  #   "value":"[OX2H0]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Acid/base",
                                                                  #   "value":"[H1,H0-]",
                                                                  #  # "description":"Works for any atom if base form has no Hs & acid has only one."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Any carbon aromatic or non-aromatic",
                                                                  #   "value":"[#6] or [c,C]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"SMILES wildcard",
                                                                  #   "value":"[#0]",
                                                                  #  # "description":"This SMARTS hits the SMILES *"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Factoring",
                                                                  #   "value":"[OX2,OX1-][OX2,OX1-] or [O;X2,X1-][O;X2,X1-]",
                                                                  #  # "description":"Factor out common atomic expressions in the recursive SMARTS. May improve human readability."
                                                                  # },
                                                                  # {
                                                                  #   "label":"High-precidence \"and\"",
                                                                  #   "value":"[N&X4&+,N&X3&+0] or [NX4+,NX3+0]",
                                                                  #  # "description":"High-precidence \"and\" (&) is the default logical operator. \"Or\" (,) is higher precidence than & and low-precidence \"and\" (;) is lower precidence than &."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Any atom w/ at-least 1 H",
                                                                  #   "value":"[*!H0,#1]",
                                                                  #  # "description":"In SMILES and SMARTS, Hydrogen is not considered an atom (unless it is specified as an isotope). The hydrogen count is instead consi dered a property of an atom. This SMARTS provides a way to effectively hit Hs themselves."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Hs on Carbons",
                                                                  #   "value":"[#6!H0,#1]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Atoms w/ 1 H",
                                                                  #   "value":"[H,#1]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Acid",
                                                                  #   "value":"[!H0;F,Cl,Br,I,N+,$([OH]-*=[!#6]),+]",
                                                                  #  # "description":"Proton donor"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Carboxylic acid",
                                                                  #   "value":"[CX3](=O)[OX2H1]",
                                                                  #  # "description":"(-oic acid, COOH)"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Carboxylic acid or conjugate base.",
                                                                  #   "value":"[CX3](=O)[OX1H0-,OX2H1]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Hydroxyl_acidic",
                                                                  #   "value":"[$([OH]-*=[!#6])]",
                                                                  #  # "description":"An acidic hydroxyl is a hydroxyl bonded to an atom which is multiply bonded to a hetero atom, this includes carboxylic, sulphur, pho sphorous, halogen and nitrogen oxyacids"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Phosphoric_Acid",
                                                                  #   "value":"[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
                                                                  #  # "description":"Hits both forms. Hits orthophosphoric acid and polyphosphoric acid anhydrides. Doesn't hit monophosphoric acid anhydride esters (in cluding acidic mono- & di- esters) but will hit some polyphosphoric acid anhydride esters (mono- esters on pyrophosphoric acid and longe r, di- esters on linear triphosphoric acid and longer). Hits acid and conjugate base."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Sulfonic Acid. High specificity.",
                                                                  #   "value":"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
                                                                  #  # "description":"Only hits carbo- sulfonic acids (Won't hit herteroatom-substituted molecules). Hits acid and conjugate base. Hits Both Depiction Fo rms. Hits Arene sulfonic acids."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Acyl Halide",
                                                                  #   "value":"[CX3](=[OX1])[F,Cl,Br,I]",
                                                                  #  # "description":"(acid halide, -oyl halide) <dl>"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Anionic divalent Nitrogen",
                                                                  #   "value":"[NX2-]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Oxenium Oxygen",
                                                                  #   "value":"[OX2H+]=*",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Oxonium Oxygen",
                                                                  #   "value":"[OX3H2+]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"Carbocation",
                                                                  #   "value":"[#6+]",
                                                                  #  # "description":""
                                                                  # },
                                                                  # {
                                                                  #   "label":"sp2 cationic carbon.",
                                                                  #   "value":"[$([cX2+](:*):*)]",
                                                                  #  # "description":"Aromatic cationic sp2 carbon with a free electron in a non-bonding sp2 hybrid orbital"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Azide ion.",
                                                                  #   "value":"[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]",
                                                                  #  # "description":"Hits N in azide ion"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Zwitterion High Specificity",
                                                                  #   "value":"[+1]~*~*~[-1]",
                                                                  #  # "description":"+1 charged atom separated by any 3 bonds from a -1 charged atom."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Zwitterion Low Specificity, Crude",
                                                                  #   "value":"[$([!-0!-1!-2!-3!-4]~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~*~*~*~[!+0!+1!+2!+3!+4])]",
                                                                  #  # "description":"Variously charged moieties separated by up to ten bonds."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Zwitterion Low Specificity",
                                                                  #   "value":"([!-0!-1!-2!-3!-4].[!+0!+1!+2!+3!+4])",
                                                                  #  # "description":"Variously charged moieties that are within the same molecule but not-necessarily connected. Uses component-level grouping. <dl>"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Hydrogen-bond acceptor",
                                                                  #   "value":"[#6,#7;R0]=[#8]",
                                                                  #  # "description":"Only hits carbonyl and nitroso. Matches a 2-atom pattern consisting of a carbon or nitrogen not in a ring, double bonded to an oxyge n."
                                                                  # },
                                                                  {
                                                                    "label":"Hydrogen-bond acceptor",
                                                                    "value":"[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]",
                                                                   # "description":"A H-bond acceptor is a heteroatom with no positive charge, note that negatively charged oxygen or sulphur are included. Excluded are halogens, including F, heteroaromatic oxygen, sulphur and pyrrole N. Higher oxidation levels of N,P,S are excluded. Note P(III) is currentl y included. Zeneca's work would imply that (O=S=O) shoud also be excluded."
                                                                  },
                                                                  {
                                                                    "label":"Hydrogen-bond donor.",
                                                                    "value":"[!$([#6,H0,-,-2,-3])]",
                                                                   # "description":"A H-bond donor is a non-negatively charged heteroatom with at least one H"
                                                                  },
                                                                  # {
                                                                  #   "label":"Hydrogen-bond donor.",
                                                                  #   "value":"[!H0;#7,#8,#9]",
                                                                  #  # "description":"Must have an N-H bond, an O-H bond, or a F-H bond"
                                                                  # },
                                                                  # {
                                                                  #   "label":"Possible intramolecular H-bond",
                                                                  #   "value":"[O,N;!H0]-*~*-*=[$([C,N;R0]=O)]",
                                                                  #  # "description":"Note that the overall SMARTS consists of five atoms. The fifth atom is defined by a \"recursive SMARTS\", where \"$()\" encloses a valid nested SMARTS and acts syntactically like an atom-primitive in the overall SMARTS. Multiple nesting is allowed."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Carbon Free-Radical",
                                                                  #   "value":"[#6;X3v3+0]",
                                                                  #  # "description":"Hits a neutral carbon with three single bonds."
                                                                  # },
                                                                  # {
                                                                  #   "label":"Nitrogen Free-Radical",
                                                                  #   "value":"[#7;X2v4+0]",
                                                                  #  # "description":"Hits a neutral nitrogen with two single bonds or with a single and a triple bond."
                                                                  # }
                                                                ], key=lambda k: k['label']) 
                                                },

                                                "dateStart": {
                                                  "title": "Added From",
                                                  "type": "string",
                                                  "format": "date",
                                                  "style": {
                                                    "margin-right": "30px;"
                                                  }, 
                                                }, 

                                                "dateEnd": {
                                                  "title": "Added Until",
                                                  "type": "string",
                                                  "format": "date",
                                                },

                                                "smiles": {
                                                  "title": "SMILES or SMARTS",
                                                  "type": "string",
                                                },

                                                "substruc": {
                                                  "title": "Structural search type",
                                                  "type": "string",
                                                  "enum": [
                                                    "with_substructure",
                                                    "flexmatch"
                                                  ],
                                                  "default": "with_substructure",
                                                },
                                                "search_custom_fields__kv_any": { 
                                                  "type": "array", 
                                                  "format" : "uiselect",
                                                  "items" :[],
                                                  "placeholder": "Tagged fields",
                                                  "title": "Any of the following custom field values:",
                                                }

                                    
                                }
                            }
                        }





    def alter_list_data_to_serialize(self, request, bundle):
        '''Here we append a list of tags to the data of the GET request if the
        search fields are required'''
        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user, request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle['user'] = userbundle.data
        editor_projects = self._meta.authorization.editor_projects(request)
        for bun in bundle["objects"]:
            bun.data["editor"] = bun.obj.id in editor_projects

        if request.GET.get("schemaform", None):
            searchfields = set([])
            searchfield_items = []

            for bun in bundle["objects"]:
                schemaform = self.get_schema_form(bun.obj.custom_field_config, 
                    bun.obj.project_key,
                    searchfield_items=searchfield_items, 
                    searchfields=searchfields,)
                bun.data["schemaform"] = schemaform
                bun.data["editor"] = bun.obj.id in editor_projects

            bundle["searchform"] = self.get_searchform(bundle,searchfield_items)
        return bundle




    def get_schema_form(self, custom_field_config, project_key, searchfield_items=[], searchfields=set([]),):
        fields = []

        for f in custom_field_config.pinned_custom_field.all():
            d = self.get_field_values(f, project_key)
            fields.append(d)
            for item in d[4]:
                if item["value"] not in searchfields:
                    searchfields.add(item["value"] )
                    searchfield_items.append(item)
        schemaform = {
                "schema" :{
                            "type" : "object",
                            "properties"   :  dict((field[0],field[1]) for field in fields),
                            "required" : []
                },
                "form" : [field[0] if not field[3] else field[3] for field in fields ]
            }
        return schemaform



    def get_field_values(self,  obj, projectKey):
        data =  copy.deepcopy(obj.FIELD_TYPE_CHOICES[obj.field_type]["data"])

        data["title"] = obj.name
        data["placeholder"] = obj.description
        data["friendly_field_type"] = obj.FIELD_TYPE_CHOICES[obj.field_type]["name"]


        form = {}
        form["field_type"] = obj.field_type
        form["position"] = obj.position
        form["key"] = obj.name
        form["title"] = obj.name
        form["placeholder"] = obj.description
        form["allowed_values"] = obj.allowed_values
        form["part_of_blinded_key"] = obj.part_of_blinded_key
        searchitems = []
        if obj.UISELECT in data.get("format", ""):
            allowed_items = obj.get_allowed_items(projectKey) 
            data["items"] = allowed_items[0]
            searchitems = allowed_items[1]
            #if we have a uiselect field with no description, make the placeholder say "Choose..."
            #if obj.description == None:
            form["placeholder"] = "Choose..."
            form["help"] = obj.description
            # form["helpdirectivename"] = "info-box"
            # form["helpdirectiveparams"] = "freetext='%s'" % (obj.description)
            # form["helpDirectiveClasses"] = "pull-right info-box"
            # #form["title"] = "%s<info-box freetext='%s'></info-box>" % (obj.name, obj.description)
        else:
            allowed_items = obj.get_allowed_items(projectKey)
            searchitems = allowed_items[1]
        


        maxdate = time.strftime("%Y-%m-%d")
        if data.get("format", False) == obj.DATE:
            form.update( {
                "minDate": "2000-01-01",
                "maxDate": maxdate,
                'type': 'datepicker',
                "format": "yyyy-mm-dd",
                'pickadate': {
                  'selectYears': True, 
                  'selectMonths': True,
                },
            })

        else:
            for item in ["options"]:
                stuff = data.pop(item, None)
                if stuff:
                    form[item] = stuff
        return (obj.name, data, obj.required, form, searchitems)

    def reindex_elasticsearch(self, request, **kwargs):
        bundle = self.build_bundle(request=request)
        #reindex compound data
        #batches = CBHCompoundBatch.objects.all()
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        es_ready_updates = [es_serializer.to_es_ready_data(proj, 
            options={"underscorize": True}) for proj in self._meta.queryset]
        index_name='chemreg_projects_index'
        elasticsearch_client.create_temporary_index(es_ready_updates, request, index_name)

        return HttpResponse(content='[]', content_type=build_content_type(desired_format) )



        
