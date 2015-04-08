from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from django.conf import settings
from django.http import HttpResponse

from tastypie.resources import ModelResource, Resource
from itertools import chain


from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch, Project, CustomFieldConfig
from cbh_chembl_ws_extension.authorization import ProjectAuthorization, ProjectListAuthorization
from tastypie.serializers import Serializer
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
import json


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
        return super(CustomFieldConfigResource, self).get_object_list(request).prefetch_related("project")

    def alter_list_data_to_serialize(self, request, bundle):
        for obj in bundle["objects"]:
            if obj.data.get("schemaform", None):
                obj.data["schemaform"] = json.loads(obj.data["schemaform"])
        return bundle
