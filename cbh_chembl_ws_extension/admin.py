from django.contrib import admin
from cbh_chembl_model_extension.models import Project, PinnedCustomField, CustomFieldConfig

from django.contrib.admin import ModelAdmin


from django.forms.widgets import HiddenInput


class GrappelliSortableHiddenMixin(object):
    """
    Mixin which hides the sortable field with Stacked and Tabular inlines.
    This mixin must precede admin.TabularInline or admin.StackedInline.
    """
    sortable_field_name = "position"

    def formfield_for_dbfield(self, db_field, **kwargs):
        if db_field.name == self.sortable_field_name:
            kwargs["widget"] = HiddenInput()
        return super(GrappelliSortableHiddenMixin, self).formfield_for_dbfield(db_field, **kwargs)


class PinnedCustomFieldInline( admin.TabularInline): #GrappelliSortableHiddenMixin
    model = PinnedCustomField
    #sortable_field_name = "position"
    extra = 0
    

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
