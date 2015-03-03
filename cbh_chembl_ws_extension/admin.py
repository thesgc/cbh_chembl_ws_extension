from django.contrib import admin
from cbh_chembl_model_extension.models import Project, PinnedCustomField, CustomFieldConfig

from django.contrib.admin import ModelAdmin

class PinnedCustomFieldInline(admin.TabularInline):
    model = PinnedCustomField


class CustomFieldConfigAdmin(ModelAdmin):

    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created' 
    inlines = [PinnedCustomFieldInline,]

    def get_changeform_initial_data(self, request):
        initial = super().get_changeform_initial_data(request)
        initial['created_by'] = request.user
        return initial


class ProjectAdmin(ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'

    def get_changeform_initial_data(self, request):
        initial = super().get_changeform_initial_data(request)
        initial['created_by'] = request.user
        return initial

admin.site.register(CustomFieldConfig, CustomFieldConfigAdmin)
admin.site.register(Project, ProjectAdmin)
