#!/usr/bin/python
# -*- coding: utf-8 -*-
from tastypie.resources import ModelResource, ALL_WITH_RELATIONS
from django.conf import settings
from django.conf.urls import url
from django.http import HttpResponse
from tastypie import fields
from cbh_core_model.models import Project
from cbh_core_ws.resources import UserResource
from cbh_core_ws.authorization import ProjectListAuthorization
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
import json
import copy
import time
from django.core.urlresolvers import reverse
from cbh_core_ws.serializers import CustomFieldsSerializer
from django.db.models import Prefetch
from cbh_core_ws.cache import CachedResource
from cbh_core_ws.resources import ProjectTypeResource, \
    CustomFieldConfigResource
from django.contrib.auth.models import User

def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """

    if 'charset' in format:
        return format
    return '%s; charset=%s' % (format, encoding)


class ChemregProjectResource(CachedResource, ModelResource):

    project_type = fields.ForeignKey(
        ProjectTypeResource, 'project_type', blank=False, null=False, full=True)
    custom_field_config = fields.ForeignKey(CustomFieldConfigResource,
                                            'custom_field_config', blank=False, null=False, full=True)
    valid_cache_get_keys = ['format', 'limit', 'project_key',
                            'schemaform']

    class Meta:

        queryset = Project.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get']
        resource_name = 'cbh_projects'
        authorization = ProjectListAuthorization()
        include_resource_uri = False
        default_format = 'application/json'

        # serializer = Serializer()

        serializer = CustomFieldsSerializer()
        filtering = {'project_key': ALL_WITH_RELATIONS}

    def get_object_list(self, request):
        return super(ChemregProjectResource,
                     self).get_object_list(request).prefetch_related(Prefetch('project_type'
                                                                              )).order_by('-modified')

    def prepend_urls(self):
        return [url(r"^(?P<resource_name>%s)/custom_fields/$"
                    % self._meta.resource_name,
                    self.wrap_view('get_custom_fields'),
                    name='get_custom_fields')]

    def get_custom_fields(self, request):
        return super(ChemregProjectResource,
                     self).get_object_list(request).prefetch_related(Prefetch('custom_field_config'
                                                                              ))

    def get_searchform(self, bundle, searchfield_items):
        '''Note that the form here is expected to have the UOx id as the first item'''
        ur = UserResource()
        uri = ur.get_resource_uri()
        return {
            'cf_form': [{
                'htmlClass': 'col-sm-10',
                'key': 'search_custom_fields__kv_any',
                'disableSuccessState': True,
                'feedback': False,
                'options': {'refreshDelay': 0,
                            'async': {'url': reverse('api_get_elasticsearch_autocomplete',
                                                     kwargs={'resource_name': 'cbh_compound_batches',
                                                             'api_name': settings.WEBSERVICES_NAME})}},
            }],
            'cf_schema': {'required': [], 'type': 'object',
                          'properties': {'search_custom_fields__kv_any': {
                              'type': 'array',
                              'format': 'uiselect',
                              'items': [],
                              'placeholder': 'Filter project data',
                              'title': 'Project data values:',
                          }}},
            'form': [
                {
                    'key': 'related_molregno__chembl__chembl_id__in',
                    'title': '%s ID' % settings.ID_PREFIX,
                    'placeholder': 'Search multiple IDs',
                    'feedback': False,
                        'htmlClass': 'col-md-4 col-xs-6',

                    'options': {'refreshDelay': 0,
                                'async': {'url': reverse('api_get_elasticsearch_ids',
                                                         kwargs={'resource_name': 'cbh_compound_batches',
                                                                 'api_name': settings.WEBSERVICES_NAME})}},
                },
                
                {
                     'key': 'creator',
                 'htmlClass': 'col-md-4 col-xs-6',
                    'placeholder': 'Select users to search',
                    'feedback': False,


                },
                {
                    'key': 'project__project_key__in',
                    'placeholder': 'Select projects to search',
                     'htmlClass': 'col-md-4 col-xs-6',
                    'feedback': False,
                    'description': 'Search for projects in order to limit the choice of fields on show.',

                                    'disableSuccessState': True,

                    'validationMessage': {'default': 'Please select a project if you wish to edit data.'}
                },
                {
                    'key': 'multiple_batch_id',
                    'htmlClass': 'col-md-4 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                },
                
                

                {
                    'key': 'dateStart',
                    'type': 'datepicker',
                    'minDate': '2004-01-01',
                    'htmlClass': 'col-md-4 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                    'pickadate': {
                        'selectYears': True,
                        'selectMonths': True
                        },
                },
                {
                    'key': 'dateEnd',
                    'type': 'datepicker',
                    'minDate': '2004-01-01',
                        'htmlClass': 'col-md-4 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                    'pickadate': {'selectYears': True,
                                  'selectMonths': True},
                },
                
                {
                    'htmlClass': 'col-md-4 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                    'key': 'functional_group',

                },
                {
                    'key': 'smiles',
                    'placeholder': 'Search SMILES or SMARTS string',
                    'append': 'today',
                    'feedback': False,
                    'htmlClass': 'col-md-4 col-xs-6',
                    'disableSuccessState': True,
                },
                {
                    'key': 'substruc',
                    'style': {'selected': 'btn-success',
                              'unselected': 'btn-default'},
                        'htmlClass': 'col-md-4 col-xs-6',
                    'type': 'radiobuttons',
                    'disableSuccessState': True,
                    'feedback': False,
                    'titleMap': [{'value': 'with_substructure',
                                  'name': 'Substructure'},
                                 {'value': 'flexmatch',
                                  'name': 'Exact Match'}],
                },
                {
                    'htmlClass': 'col-md-4 col-xs-6',
                    'key': 'search_custom_fields__kv_any',
                    'disableSuccessState': True,
                    'help': 'Searching using this filter will bring back results that match an OR pattern within the same data category, with AND across data categories, i.e. results which contain this item within category a OR that item within category a AND that item within category b.',
                    'feedback': False,
                    'options': {'refreshDelay': 0,
                                'async': {'url': reverse('api_get_elasticsearch_autocomplete',
                                                         kwargs={'resource_name': 'cbh_compound_batches',
                                                                 'api_name': settings.WEBSERVICES_NAME})},
                                "tagging": "tagFunction",
                                "taggingLabel": "(in any field)",
                                "taggingTokens": "",
                                },
                },

            ],
            'schema': {'required': [], 'type': 'object', 'properties': {
                'related_molregno__chembl__chembl_id__in': {
                    'type': 'array',
                    'format': 'uiselect',
                    
                },

                'creator': {
                    'type': 'array',
                    'format': 'uiselect',
                       'title': 'Compound batch created by',
                        'type': 'array',
                        'format': 'uiselect',
                        'htmlClass': 'col-md-4 col-xs-6',
                        'placeholder': 'Search user who created the batch',
                        'options': {'searchDescriptions': False},
                        'items':  sorted([
                            {'label': user.first_name + " " + user.last_name + ""+ user.username + "", "value" : uri + '/' + str(user.id) } 
                            for user in User.objects.exclude(pk=-1)
                        ], key=lambda k: k['label'])
                },
                
                'multiple_batch_id': {'title': 'Upload ID',
                                      'type': 'string'},
                'project__project_key__in': {
                    'title': 'Project',
                    'type': 'array',
                    'format': 'uiselect',
                    'items': [{'label': p.obj.name,
                               'value': p.obj.project_key} for p in
                              bundle['objects']],
                },
                'functional_group': {
                    'title': 'Functional Group',
                    'type': 'string',
                    'format': 'uiselect',
                    'placeholder': 'Search chemical groups',
                    'options': {'searchDescriptions': False},
                    'default': '',
                    'copyValueTo': 'smiles',
                    'items': [{'label': 'None', 'value': ''}] + sorted([
                        {'label': 'Alkyl Carbon', 'value': '[CX4]'},
                        {'label': 'Allenic Carbon',
                         'value': '[$([CX2](=C)=C)]'},
                        {'label': 'Vinylic Carbon',
                         'value': '[$([CX3]=[CX3])]'},
                        {'label': 'Acetylenic Carbon',
                         'value': '[$([CX2]#C)]'},
                        {'label': 'Arene', 'value': 'c'},
                        {'label': 'Carbonyl group. Low specificity',
                         'value': '[CX3]=[OX1]'},
                        {'label': 'Carbonyl group',
                         'value': '[$([CX3]=[OX1]),$([CX3+]-[OX1-])]'},
                        {'label': 'Carbonyl with Carbon',
                         'value': '[CX3](=[OX1])C'},
                        {'label': 'Carbonyl with Nitrogen.',
                         'value': '[OX1]=CN'},
                        {'label': 'Carbonyl with Oxygen.',
                         'value': '[CX3](=[OX1])O'},
                        {'label': 'Acyl Halide',
                         'value': '[CX3](=[OX1])[F,Cl,Br,I]'},
                        {'label': 'Aldehyde', 'value': '[CX3H1](=O)[#6]'
                         },
                        {'label': 'Anhydride',
                         'value': '[CX3](=[OX1])[OX2][CX3](=[OX1])'},
                        {'label': 'Amide',
                         'value': '[NX3][CX3](=[OX1])[#6]'},
                        {'label': 'Amidinium',
                         'value': '[NX3][CX3]=[NX3+]'},
                        {'label': 'Carbamate.',
                         'value': '[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]'},
                        {'label': 'Carbamic ester',
                         'value': '[NX3][CX3](=[OX1])[OX2H0]'},
                        {'label': 'Carbamic acid.',
                         'value': '[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]'
                         },
                        {'label': 'Carboxylate Ion.',
                         'value': '[CX3](=O)[O-]'},
                        {'label': 'Carbonic Acid or Carbonic Ester',
                         'value': '[CX3](=[OX1])(O)O'},
                        {'label': 'Carbonic Acid or Carbonic Acid-Ester', 'value': '[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]'
                         },
                        {'label': 'Carbonic Ester (carbonic acid diester)',
                         'value': 'C[OX2][CX3](=[OX1])[OX2]C'},
                        {'label': 'Carboxylic acid',
                         'value': '[CX3](=O)[OX2H1]'},
                        {'label': 'Carboxylic acid or conjugate base.',
                         'value': '[CX3](=O)[OX1H0-,OX2H1]'},
                        {'label': 'Cyanamide',
                         'value': '[NX3][CX2]#[NX1]'},
                        {'label': 'Ester Also hits anhydrides',
                         'value': '[#6][CX3](=O)[OX2H0][#6]'},
                        {'label': 'Ketone', 'value': '[#6][CX3](=O)[#6]'
                         },
                        {'label': 'Ether', 'value': '[OD2]([#6])[#6]'},
                        {'label': 'Hydrogen Atom', 'value': '[H]'},
                        {'label': 'Not a Hydrogen Atom',
                         'value': '[!#1]'},
                        {'label': 'Proton', 'value': '[H+]'},
                        {'label': 'Mono-Hydrogenated Cation',
                         'value': '[+H]'},
                        {'label': 'Not Mono-Hydrogenated',
                         'value': '[!H] or [!H1]'},
                        {'label': 'Primary or secondary amine, not amide.',
                            'value': '[NX3;H2,H1;!$(NC=O)]'},
                        {'label': 'Enamine', 'value': '[NX3][CX3]=[CX3]'
                         },
                        {'label': 'Primary amine, not amide.',
                         'value': "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6] Not amide (C not double bonded to a hetero-atom), not ammonium ion (N must be 3-connected), not ammonia (N's H-count can't be 3), not cyanamide (C not triple bonded to a hetero-atom)"
                         },
                        {'label': 'Two primary or secondary amines',
                         'value': '[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]'
                         },
                        {'label': 'Enamine or Aniline Nitrogen',
                         'value': '[NX3][$(C=C),$(cc)]'},
                        {'label': 'Azide group.',
                         'value': '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]'
                         },
                        {'label': 'Azide ion.',
                         'value': '[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]'
                         },
                        {'label': 'Nitrogen.', 'value': '[#7]'},
                        {'label': 'Azo Nitrogen. Low specificity.',
                         'value': '[NX2]=N'},
                        {'label': 'Azo Nitrogen.diazene',
                         'value': '[NX2]=[NX2]'},
                        {'label': 'Azoxy Nitrogen.',
                         'value': '[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]'
                         },
                        {'label': 'Diazo Nitrogen',
                         'value': '[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]'
                         },
                        {'label': 'Azole.',
                         'value': '[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]'
                         },
                        {'label': 'Hydrazine H2NNH2',
                         'value': '[NX3][NX3]'},
                        {'label': 'Hydrazone C=NNH2',
                         'value': '[NX3][NX2]=[*]'},
                        {'label': 'Substituted imine',
                         'value': '[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]'
                         },
                        {'label': 'Substituted or un-substituted imine',
                         'value': '[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]'
                         },
                        {'label': 'Iminium', 'value': '[NX3+]=[CX3]'},
                        {'label': 'Unsubstituted dicarboximide',
                         'value': '[CX3](=[OX1])[NX3H][CX3](=[OX1])'},
                        {'label': 'Substituted dicarboximide',
                         'value': '[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])'
                         },
                        {'label': 'Dicarboxdiimide',
                         'value': '[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])'
                         },
                        {'label': 'Nitrate group',
                         'value': '[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]'
                         },
                        {'label': 'Nitrate Anion',
                         'value': '[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]'
                         },
                        {'label': 'Nitrile', 'value': '[NX1]#[CX2]'},
                        {'label': 'Isonitrile', 'value': '[CX1-]#[NX2+]'
                         },
                        {'label': 'Nitro group.',
                         'value': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8] Hits both forms.'
                         },
                        {'label': 'Two Nitro groups',
                         'value': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8].[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]'
                         },
                        {'label': 'Nitroso-group',
                         'value': '[NX2]=[OX1]'},
                        {'label': 'N-Oxide',
                         'value': '[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]'
                         },
                        {'label': 'Hydroxyl', 'value': '[OX2H]'},
                        {'label': 'Hydroxyl in Alcohol',
                         'value': '[#6][OX2H]'},
                        {'label': 'Hydroxyl in Carboxylic Acid',
                         'value': '[OX2H][CX3]=[OX1]'},
                        {'label': 'Hydroxyl in H-O-P-',
                         'value': '[OX2H]P'},
                        {'label': 'Enol', 'value': '[OX2H][#6X3]=[#6]'
                         },
                        {'label': 'Phenol', 'value': '[OX2H][cX3]:[c]'
                         },
                        {'label': 'Enol or Phenol',
                         'value': '[OX2H][$(C=C),$(cc)]'},
                        {'label': 'Hydroxyl_acidic',
                         'value': '[$([OH]-*=[!#6])]'},
                        {'label': 'Peroxide groups.',
                         'value': '[OX2,OX1-][OX2,OX1-]'},
                        {'label': 'Phosphoric_acid groups.',
                         'value': '[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]'
                         },
                        {'label': 'Phosphoric_ester groups.',
                         'value': '[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]'
                         },
                        {'label': 'Carbo-Thiocarboxylate',
                         'value': '[S-][CX3](=S)[#6]'},
                        {'label': 'Carbo-Thioester',
                         'value': 'S([#6])[CX3](=O)[#6]'},
                        {'label': 'Thio analog of carbonyl',
                         'value': '[#6X3](=[SX1])([!N])[!N]'},
                        {'label': 'Thiol, Sulfide or Disulfide Sulfur',
                         'value': '[SX2]'},
                        {'label': 'Thiol', 'value': '[#16X2H]'},
                        {'label': 'Sulfur with at-least one hydrogen.',
                         'value': '[#16!H0]'},
                        {'label': 'Thioamide',
                         'value': '[NX3][CX3]=[SX1]'},
                        {'label': 'Sulfide', 'value': '[#16X2H0]'},
                        {'label': 'Mono-sulfide',
                         'value': '[#16X2H0][!#16]'},
                        {'label': 'Di-sulfide',
                         'value': '[#16X2H0][#16X2H0]'},
                        {'label': 'Two Sulfides',
                         'value': '[#16X2H0][!#16].[#16X2H0][!#16]'},
                        {'label': 'Sulfinate',
                         'value': '[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]'
                         },
                        {'label': 'Sulfinic Acid',
                         'value': '[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]'
                         },
                        {'label': 'Sulfone. Low specificity.',
                         'value': '[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]'
                         },
                        {'label': 'Sulfone. High specificity.',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]'
                         },
                        {'label': 'Sulfonic acid. High specificity.',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]'
                         },
                        {'label': 'Sulfonate',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]'
                         },
                        {'label': 'Sulfonamide.',
                         'value': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]'
                         },
                        {'label': 'Carbo-azosulfone',
                         'value': '[SX4](C)(C)(=O)=N'},
                        {'label': 'Sulfonamide',
                         'value': '[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]'
                         },
                        {'label': 'Sulfoxide Low specificity.',
                         'value': '[$([#16X3]=[OX1]),$([#16X3+][OX1-])]'
                         },
                        {'label': 'Sulfoxide High specificity',
                         'value': '[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]'
                         },
                        {'label': 'Sulfate',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]'
                         },
                        {'label': 'Sulfuric acid ester (sulfate ester) Low specificity.',
                         'value': '[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]'
                         },
                        {'label': 'Sulfuric Acid Diester.',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]'
                         },
                        {'label': 'Sulfamate.',
                         'value': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]'
                         },
                        {'label': 'Sulfamic Acid.',
                         'value': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]'
                         },
                        {'label': 'Sulfenic acid.',
                         'value': '[#16X2][OX2H,OX1H0-]'},
                        {'label': 'Sulfenate.',
                         'value': '[#16X2][OX2H0]'},
                        {'label': 'Any carbon attached to any halogen',
                         'value': '[#6][F,Cl,Br,I]'},
                        {'label': 'Halogen', 'value': '[F,Cl,Br,I]'},
                        {'label': 'Sulfide', 'value': '[#16X2H0]'},
                        {'label': 'Mono-sulfide',
                         'value': '[#16X2H0][!#16]'},
                        {'label': 'Di-sulfide',
                         'value': '[#16X2H0][#16X2H0]'},
                        {'label': 'Hydrogen-bond acceptor',
                         'value': '[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'
                         },
                        {'label': 'Hydrogen-bond donor.',
                         'value': '[!$([#6,H0,-,-2,-3])]'},
                    ], key=lambda k: k['label']),
                },
                'dateStart': {
                    'title': 'Added after',
                    'type': 'string',
                    'format': 'date',
                    'style': {'margin-right': '30px;'},
                },
                'dateEnd': {'title': 'Added before', 'type': 'string',
                            'format': 'date'},
                'smiles': {'title': 'SMILES or SMARTS',
                           'type': 'string'},
                'substruc': {
                    'title': 'Structural search type',
                    'type': 'string',
                    'enum': ['with_substructure', 'flexmatch'],
                    'default': 'with_substructure',
                },
                'search_custom_fields__kv_any': {
                    'type': 'array',
                    'format': 'uiselect',
                    'items': [],
                    'placeholder': 'Choose column and value...',
                    'title': 'Filter by project data values:',
                },

            }},
        }

    def alter_list_data_to_serialize(self, request, bundle):
        '''Here we append a list of tags to the data of the GET request if the
        search fields are required'''

        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user,
                                          request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle['user'] = userbundle.data

        editor_projects = \
            self._meta.authorization.editor_projects(request)
        for bun in bundle['objects']:
            bun.data['editor'] = bun.obj.id in editor_projects

        if request.GET.get('schemaform', None):
            searchfields = set([])
            searchfield_items = []

            for bun in bundle['objects']:
                schemaform = \
                    self.get_schema_form(bun.obj.custom_field_config,
                                         bun.obj.project_key,
                                         searchfield_items=searchfield_items,
                                         searchfields=searchfields)
                bun.data['schemaform'] = schemaform
                bun.data['editor'] = bun.obj.id in editor_projects

            bundle['searchform'] = self.get_searchform(bundle,
                                                       searchfield_items)

        if self.determine_format(request) \
                == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' \
                or request.GET.get('format') == 'xls':

            cfr_string = \
                self.get_object_list(request).filter(id=request.GET.get('project_key'
                                                                        ))[0].custom_field_config.schemaform
            cfr_json = json.loads(cfr_string)
            bundle['custom_field_config'] = cfr_json['form']

        return bundle

    def get_schema_form(
        self,
        custom_field_config,
        project_key,
        searchfield_items=[],
        searchfields=set([]),
    ):
        fields = []

        for f in custom_field_config.pinned_custom_field.all():
            d = self.get_field_values(f, project_key)
            fields.append(d)
            for item in d[4]:
                if item['value'] not in searchfields:
                    searchfields.add(item['value'])
                    searchfield_items.append(item)
        schemaform = {'schema': {'type': 'object',
                                 'properties': dict((field[0], field[1])
                                                    for field in fields), 'required': []},
                      'form': [(field[0] if not field[3] else field[3])
                               for field in fields]}
        return schemaform

    def get_field_values(self, obj, projectKey):
        data = \
            copy.deepcopy(obj.FIELD_TYPE_CHOICES[obj.field_type]['data'
                                                                 ])

        data['title'] = obj.name
        data['placeholder'] = obj.description
        data['friendly_field_type'] = \
            obj.FIELD_TYPE_CHOICES[obj.field_type]['name']

        form = {}
        form['field_type'] = obj.field_type
        form['position'] = obj.position
        form['key'] = obj.name
        form['title'] = obj.name
        form['placeholder'] = obj.description
        form['allowed_values'] = obj.allowed_values
        form['part_of_blinded_key'] = obj.part_of_blinded_key
        searchitems = []
        if obj.UISELECT in data.get('format', ''):
            allowed_items = obj.get_allowed_items(projectKey)
            data['items'] = allowed_items[0]
            searchitems = allowed_items[1]

            # if we have a uiselect field with no description, make the placeholder say "Choose..."
            # if obj.description == None:

            form['placeholder'] = 'Choose...'
            form['help'] = obj.description
        else:

            # form["helpdirectivename"] = "info-box"
            # form["helpdirectiveparams"] = "freetext='%s'" % (obj.description)
            # form["helpDirectiveClasses"] = "pull-right info-box"
            # #form["title"] = "%s<info-box freetext='%s'></info-box>" % (obj.name, obj.description)

            allowed_items = obj.get_allowed_items(projectKey)
            searchitems = allowed_items[1]

        maxdate = time.strftime('%Y-%m-%d')
        if data.get('format', False) == obj.DATE:
            form.update({
                'minDate': '2000-01-01',
                'maxDate': maxdate,
                'type': 'datepicker',
                'format': 'yyyy-mm-dd',
                'pickadate': {'selectYears': True,
                              'selectMonths': True},
            })
        else:

            for item in ['options']:
                stuff = data.pop(item, None)
                if stuff:
                    form[item] = stuff
        return (obj.name, data, obj.required, form, searchitems)

    def create_response(
        self,
        request,
        data,
        response_class=HttpResponse,
        **response_kwargs
    ):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """

        desired_format = self.determine_format(request)
        serialized = self.serialize(request, data, desired_format)
        rc = response_class(content=serialized,
                            content_type=build_content_type(desired_format),
                            **response_kwargs)

        if desired_format \
                == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
            rc['Content-Disposition'] = \
                'attachment; filename=project_data_explanation.xlsx'
        return rc
