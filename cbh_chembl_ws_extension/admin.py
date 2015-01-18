from django.contrib import admin
from cbh_chembl_model_extension.models import Project

from guardian.admin import GuardedModelAdmin

class ProjectAdmin(GuardedModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'

admin.site.register(Project, ProjectAdmin)
