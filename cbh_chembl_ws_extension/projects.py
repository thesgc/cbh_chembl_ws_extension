from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from django.conf import settings
from django.http import HttpResponse

from tastypie.resources import ModelResource, Resource
from itertools import chain


from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch, Project, CustomFieldConfig, PinnedCustomField
from cbh_chembl_ws_extension.authorization import ProjectAuthorization, ProjectListAuthorization
from tastypie.serializers import Serializer
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
import json
import copy
import time


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

    def get_searchform(self, bundle,searchfield_items ):
        return { "form": [
                                {
                                  "key": "project",
                                  "type": "checkboxes",
                                  "description": "<info-box freetext='Limit your search results to items tagged with project-specific information'></info-box>",
                                  "placeholder": "Check the console",
                                  "htmlClass": "col-xs-12",
                                  "onChange": "getSearchCustomFields()",
                                  "titleMap": {
                                        p.obj.project_key : p.obj.name for p in  bundle["objects"]
                                  },
                                  
                                  "disableSuccessState": True,
                                },
                                {
                                  "key": "dateStart",
                                  "minDate": "2004-01-01",
                                  "htmlClass": "col-xs-6",
                                  "disableSuccessState": True,
                                },
                                {
                                  "key": "dateEnd",
                                  "minDate": "2004-01-01",
                                  "htmlClass": "col-xs-6",
                                  "disableSuccessState": True,
                                },
                                {
                                  "key": "smiles",
                                  "placeholder": "Search SMILES string",
                                  "append": "today",
                                  "feedback": False,
                                  "htmlClass": "col-xs-12",
                                  "disableSuccessState": True,
                                },
                                {
                                  "key": "substruc",
                                  "style": {
                                    "selected": "btn-success",
                                    "unselected": "btn-default"
                                  },
                                  "htmlClass": "col-xs-12",
                                  "type": "radiobuttons",
                                  "disableSuccessState": True,
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
                                    "type": "fieldset",
                                    "htmlClass": "col-xs-12",
                                    "items": [
                                     "search_custom_fields__kv_any",
                                ]
                                }
                            ],
                            "schema": {
                                "required": [
                                ],
                                "type": "object",
                                "properties": {
                                                "project": {
                                                  "title": "Project",
                                                  "type": "string",
                                                    "enum" :[p.obj.project_key for p in bundle["objects"]],
                                                 "default" : [p.obj.project_key for p in bundle["objects"]]

                                                },

                                                "dateStart": {
                                                  "title": "From",
                                                  "type": "string",
                                                  "format": "date",
                                                  "style": {
                                                    "margin-right": "30px;"
                                                  }, 
                                                }, 

                                                "dateEnd": {
                                                  "title": "To",
                                                  "type": "string",
                                                  "format": "date",
                                                },

                                                "smiles": {
                                                  "title": "SMILES",
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
                                                        "items" : sorted(searchfield_items),
                                                         "placeholder": "Tagged fields",
                                                         "title": "Limit by the following tags:",
                                                   
                                                }

                                    
                                }
                            }
                        }





    def alter_list_data_to_serialize(self, request, bundle):
        '''Here we append a list of tags to the data of the GET request if the
        search fields are required'''
        if request.GET.get("schemaform", None):
            searchfields = set([])
            searchfield_items = []

            for bun in bundle["objects"]:
                schemaform = self.get_schema_form(bun.obj.custom_field_config, 
                    bun.obj.project_key,
                    searchfield_items=searchfield_items, 
                    searchfields=searchfields,)
                bun.data["schemaform"] = schemaform

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


        form = {}
        form["field_type"] = obj.field_type
        form["positon"] = obj.position
        form["key"] = obj.name
        form["title"] = obj.name
        form["placeholder"] = obj.description
        form["allowed_values"] = obj.allowed_values
        searchitems = []
        if obj.UISELECT in data.get("format", ""):
            allowed_items = obj.get_allowed_items(projectKey) 
            data["items"] = allowed_items[0]
            searchitems = allowed_items[1]

          
        if data.get("format", False) == obj.DATE:
            form.update( {
                "minDate": "2000-01-01",
                "maxDate": time.strftime("%d-%m-%Y"),
                "format": "dd-mm-yyyy"

            })

        else:
            for item in ["options"]:
                stuff = data.pop(item, None)
                if stuff:
                    form[item] = stuff
        return (obj.name, data, obj.required, form, searchitems)



        
