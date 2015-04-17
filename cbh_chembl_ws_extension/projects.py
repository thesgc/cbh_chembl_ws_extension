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



class PinnedCustomFieldResource(ModelResource):
    class Meta:
        queryset = PinnedCustomField.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get']        
        resource_name = 'cbh_pinned_custom_fields'
        #authorization = ProjectAuthorization()
        include_resource_uri = False
        default_format = 'application/json'
        serializer = Serializer()
        
    def alter_detail_data_to_serialize(self, request, bundle):
        bundle = bundle.obj.get_allowed_items(request.GET.get("projectKey",""))
        return bundle

    def get_field_values(self, obj, projectKey):
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
        if data.get("format", False) == obj.UISELECT:
            data["items"] = obj.get_allowed_items(projectKey) 

          
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
        return (obj.name, data, obj.required, form)

class CustomFieldConfigResource(ModelResource):


    class Meta:
        queryset = CustomFieldConfig.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get']        
        resource_name = 'cbh_custom_field_configs'
        authorization = ProjectAuthorization()
        include_resource_uri = False
        default_format = 'application/json'
        serializer = Serializer()
        filtering = {
            
            "project__project_key": ALL_WITH_RELATIONS,
        }


    def get_object_list(self, request):
        return super(CustomFieldConfigResource, self).get_object_list(request).prefetch_related("project").select_related("pinned_custom_field")



    def alter_list_data_to_serialize(self, request, bundle):
        for bun in bundle["objects"]:
            bun.data["schemaform"] = self.get_schema_form(bun.obj, request.GET.get("projectKey", ""))
        return bundle

    def get_schema_form(self, obj, projectKey):
        fields = []
        # self.select_related()
        for f in obj.pinned_custom_field.all():
            pcfr = PinnedCustomFieldResource()
            d = pcfr.get_field_values(f, projectKey)
            fields.append(d)
        schemaform = {
                        "schema" :{
                                    "type" : "object",
                                    "properties"   :  dict((field[0],field[1]) for field in fields),
                                    "required" : []
                        },
                        "form" : [field[0] if not field[3] else field[3] for field in fields ]
                    }
        return schemaform