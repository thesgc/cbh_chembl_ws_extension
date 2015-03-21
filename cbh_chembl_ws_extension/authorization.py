from tastypie.authorization import Authorization
from tastypie.exceptions import TastypieError, Unauthorized
import logging
logger = logging.getLogger(__name__)
logger_debug = logging.getLogger(__name__)
from cbh_chembl_model_extension.models import Project





class ProjectListAuthorization(Authorization):
    """
    Uses permission checking from ``django.contrib.auth`` to map
    ``POST / PUT / DELETE / PATCH`` to their equivalent Django auth
    permissions.

    Both the list & detail variants simply check the model they're based
    on, as that's all the more granular Django's permission setup gets.
    """

    def login_checks(self, request, model_klass):

        # If it doesn't look like a model, we can't check permissions.
        # if not model_klass or not getattr(model_klass, '_meta', None):
        #     print "improper_setup_of_authorization"
        #     raise Unauthorized("improper_setup_of_authorization")
        # User must be logged in to check permissions.
        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")


    def base_checks(self, request, model_klass, data, possible_perm_levels):
        self.login_checks(request, model_klass)
        if not data.get("project_key", None):
            print "no_project_key"
            raise Unauthorized("no_project_key")

        project = data.get("project", False)
        
        has_perm = Project.objects.get_user_permission(project.id, request.user, possible_perm_levels)
        if has_perm is True:
            return True

        print "user_does_not_have_correct_permissions_for_operation"
        raise Unauthorized("user_does_not_have_correct_permissions_for_operation")
        

    def list_checks(self, request, model_klass, data, possible_perm_levels, object_list):
        self.login_checks(request, model_klass)
        new_list = []
        for obj in object_list:
            if Project.objects.get_user_permission(obj.id, request.user, possible_perm_levels) is True:
                new_list.append(obj)

        return new_list



    def read_list(self, object_list, bundle):
        return self.list_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor","viewer",], object_list)

        


















class ProjectAuthorization(Authorization):
    """
    Uses permission checking from ``django.contrib.auth`` to map
    ``POST / PUT / DELETE / PATCH`` to their equivalent Django auth
    permissions.

    Both the list & detail variants simply check the model they're based
    on, as that's all the more granular Django's permission setup gets.
    """

    def login_checks(self, request, model_klass):

        # If it doesn't look like a model, we can't check permissions.
        # if not model_klass or not getattr(model_klass, '_meta', None):
        #     print "improper_setup_of_authorization"
        #     raise Unauthorized("improper_setup_of_authorization")
        # User must be logged in to check permissions.
        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")


    def base_checks(self, request, model_klass, data, possible_perm_levels):
        self.login_checks(request, model_klass)
        if not data.get("project_key", None):
            print "no_project_key"
            raise Unauthorized("no_project_key")

        project = data.get("project", False)
        
        has_perm = Project.objects.get_user_permission(project.id, request.user, possible_perm_levels)
        if has_perm is True:
            return True

        print "user_does_not_have_correct_permissions_for_operation"
        raise Unauthorized("user_does_not_have_correct_permissions_for_operation")
        

    def list_checks(self, request, model_klass, data, possible_perm_levels, object_list):
        self.login_checks(request, model_klass)
        new_list = []
        for obj in object_list:
            if Project.objects.get_user_permission(obj.project_id, request.user, possible_perm_levels) is True:
                new_list.append(obj)

        return new_list





    #     if klass is False:
    #         return []
    #     permission = '%s.view_%s' % (klass._meta.app_label, klass._meta.module_name)
    #     for obj in object_list:        
    #         #if bundle.request.user.has_perms(permission,obj):
    #             read_list.append(obj)
    #     # GET-style methods are always allowed.
    #     if read_list:
    #         return read_list
    #     raise Unauthorized("You are not allowed to access that resource.")

    # def read_detail(self, object_list, bundle):
    #     klass = self.base_checks(bundle.request, bundle.obj.__class__, bundle.data)
    #     read_list=[]


    #     if klass is False:
    #         raise Unauthorized("You are not allowed to access that resource.")

    #     permission = '%s.view_%s' % (klass._meta.app_label, klass._meta.module_name)
    #     for obj in object_list:        
    #         #if bundle.request.user.has_perms(permission,obj):
    #             read_list.append(obj)
                
    #     if read_list:
    #         return True
    #     raise Unauthorized("You are not allowed to access that resource.")

    def create_list(self, object_list, bundle):
        bool = self.base_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor",])
        if bool is True:
            return object_list
        else:
            
            return []


    def read_detail(self, object_list, bundle):
        self.login_checks(bundle.request, bundle.obj.__class__,)
        return Project.objects.get_user_permission(bundle.obj.project_id, bundle.request.user, ["editor","viewer",]) 
                #return self.base_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor","viewer",])


    
    def create_detail(self, object_list, bundle):
        return self.base_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor",])


    def update_detail(self, object_list, bundle):
        return self.base_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor",])


    def read_list(self, object_list, bundle):
        return self.list_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor","viewer",], object_list)

        

        # if klass is False:
        #     raise Unauthorized("You are not allowed to access that resource.")

        # permission = '%s.add_%s' % (klass._meta.app_label, klass._meta.module_name)

        # for obj in object_list:        
        #     #if bundle.request.user.has_perms(permission,obj):
        #         create_list.append(obj)
                
        # if create_list:
        #     return True
        # raise Unauthorized("You are not allowed to access that resource.")

    # def update_list(self, object_list, bundle):
    #     klass = self.base_checks(bundle.request, object_list.model, bundle.data)
    #     update_list=[]

    #     if klass is False:
    #         return []

    #     permission = '%s.change_%s' % (klass._meta.app_label, klass._meta.module_name)

    #     for obj in object_list:        
    #         #if bundle.request.user.has_perms(permission,obj):
    #             update_list.append(obj)

    #     if update_list:
    #         return update_list
    #     raise Unauthorized("You are not allowed to access that resource.")

    # def update_detail(self, object_list, bundle):
    #     update_list=[]
    #     klass = self.base_checks(bundle.request, bundle.obj.__class__, bundle.data)

    #     if klass is False:
    #         raise Unauthorized("You are not allowed to access that resource.")

    #     permission = '%s.change_%s' % (klass._meta.app_label, klass._meta.module_name)

    #     for obj in object_list:        
    #         #if bundle.request.user.has_perms(permission,obj):
    #             update_list.append(obj)

    #     if update_list:
    #         return update_list
    #     raise Unauthorized("You are not allowed to access that resource.")

    # def delete_list(self, object_list, bundle):
    #     delete_list=[]
    #     klass = self.base_checks(bundle.request, object_list.model, bundle.data)

    #     if klass is False:
    #         return []

    #     permission = '%s.delete_%s' % (klass._meta.app_label, klass._meta.module_name)

    #     for obj in object_list:        
    #         #if bundle.request.user.has_perms(permission,obj):
    #             delete_list.append(obj)

    #     if delete_list:
    #         return delete_list
    #     raise Unauthorized("You are not allowed to access that resource.")

    # def delete_detail(self, object_list, bundle):
    #     delete_list=[]

    #     klass = self.base_checks(bundle.request, bundle.obj.__class__, bundle.data)

    #     if klass is False:
    #         raise Unauthorized("You are not allowed to access that resource.")

    #     permission = '%s.delete_%s' % (klass._meta.app_label, klass._meta.module_name)

    #     for obj in object_list:        
    #         #if bundle.request.user.has_perms(permission,obj):
    #             delete_list.append(obj)

    #     if delete_list:
    #         return delete_list
    #     raise Unauthorized("You are not allowed to access that resource.")
