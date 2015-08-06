from tastypie.authorization import Authorization
from tastypie.exceptions import TastypieError, Unauthorized
import logging
logger = logging.getLogger(__name__)
logger_debug = logging.getLogger(__name__)
from cbh_core_model.models import Project

















def get_all_project_ids_for_user_perms( perms, possible_perm_levels ):
    pids = []
    for perm in perms:
        prms= str(perm).split(".")
        pid = prms[0]
        if pid[0].isdigit()  and prms[1] in possible_perm_levels:
           pids.append(int(pid))
    return pids



class ProjectListAuthorization(Authorization):
    """
    Uses permission checking from ``django.contrib.auth`` to map
    ``POST / PUT / DELETE / PATCH`` to their equivalent Django auth
    permissions.

    Both the list & detail variants simply check the model they're based
    on, as that's all the more granular Django's permission setup gets.
    """
    def editor_projects(self, request, ):
        pids = get_all_project_ids_for_user_perms( request.user.get_all_permissions(), ["editor"] )
        return pids
        
    def login_checks(self, request, model_klass):

        # If it doesn't look like a model, we can't check permissions.
        # if not model_klass or not getattr(model_klass, '_meta', None):
        #     print "improper_setup_of_authorization"
        #     raise Unauthorized("improper_setup_of_authorization")
        # User must be logged in to check permissions.
        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")


    def base_checks(self, request, model_klass, data, possible_perm_levels, perms=None):
        self.login_checks(request, model_klass)
        if not data.get("project_key", None):
            print "no_project_key"
            raise Unauthorized("no_project_key")

        project = data.get("project", False)
        
        has_perm = Project.objects.get_user_permission(project.id, request.user, possible_perm_levels, perms=perms)
        if has_perm is True:
            return True

        print "user_does_not_have_correct_permissions_for_operation"
        raise Unauthorized("user_does_not_have_correct_permissions_for_operation")
        

    def list_checks(self, request, model_klass, data, possible_perm_levels, object_list):
        perms = request.user.get_all_permissions()
        pids = get_all_project_ids_for_user_perms( perms, possible_perm_levels )
        self.login_checks(request,  model_klass, )
    
        return object_list.filter(pk__in=pids)

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


    def login_checks(self, request, model_klass, perms=None):

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

        if not data.get("project__project_key", None):
            if not data.get("project_key"):
                if not data.get("projectKey"):
                    try:
                        key = data.project.project_key
                    except:

                        print "no_project_key"
                        raise Unauthorized("no_project_key")
                else:
                    key=data.get("projectKey")
            else:
                key=data.get("project_key")
        else:
            key = data.get("project__project_key")

        project = Project.objects.get(project_key=key)
        pids = get_all_project_ids_for_user_perms(request.user.get_all_permissions(), possible_perm_levels)
        if project.id in pids:
            return True
        return False
        

    def project_ids(self, request, ):
        pids = get_all_project_ids_for_user_perms( request.user.get_all_permissions(), ["editor","viewer",] )
        return pids



    def create_list(self, object_list, bundle):
        print "create"
        bool = self.base_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor",])
        if bool is True:
            return object_list
        else:
            
            return []


    def read_detail(self, object_list, bundle):
        print "readdet"

        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = get_all_project_ids_for_user_perms(bundle.request.user.get_all_permissions(), ["editor","viewer"])
        if bundle.obj.project.id in pids:
            return True
        else:
            raise Unauthorized("not authroized for project")



    def update_list(self, object_list, bundle):
        print "update"

        return []


    def create_detail(self, object_list, bundle):
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = get_all_project_ids_for_user_perms(bundle.request.user.get_all_permissions(), ["editor",])
        if bundle.data["project"].id in pids:
            return True
        else:
            raise Unauthorized("not authroized for project")
        #return self.base_checks(bundle.request, bundle.obj.__class__, bundle.data, ["editor",])


    def update_detail(self, object_list, bundle):
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = get_all_project_ids_for_user_perms(bundle.request.user.get_all_permissions(), ["editor",])
        if bundle.obj.project.id in pids:
            return True

        raise Unauthorized("not authroized for project")


    def read_list(self, object_list, bundle):
        return object_list



        
