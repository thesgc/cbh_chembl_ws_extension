# -*- coding: utf-8 -*-

from django.contrib import admin
from cbh_chembl_model_extension.models import Project, PinnedCustomField, CustomFieldConfig, SkinningConfig

from django.contrib.admin import ModelAdmin
from cbh_chembl_ws_extension.projects import ProjectResource
from cbh_chembl_ws_extension.compounds import CBHCompoundBatchResource

from django.forms.widgets import HiddenInput, TextInput
from django.db import models
import json
from solo.admin import SingletonModelAdmin

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


class PinnedCustomFieldInline( GrappelliSortableHiddenMixin, admin.TabularInline, ): #GrappelliSortableHiddenMixin
    model = PinnedCustomField
    exclude = ["field_key"]

    sortable_field_name = "position"
    formfield_overrides = {
        models.CharField: {'widget': TextInput(attrs={'size':'20'})},
    }
    extra = 3
    def get_extra (self, request, obj=None, **kwargs):
        """Dynamically sets the number of extra forms. 0 if the related object
        already exists or the extra configuration otherwise."""
        if obj:
            # Don't add any extra forms if the related object already exists.
            return 0
        return self.extra
#Make a template have to be chosen in order to create a schema and make it impossible to edit schemas once created then versioning not needed



class CustomFieldConfigAdmin(ModelAdmin):
    
    exclude= ["created_by", ]

    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created' 
    inlines = [PinnedCustomFieldInline,]

    
    def get_readonly_fields(self, request, obj=None):
        if obj: # editing an existing object
            return self.readonly_fields + ('schemaform',)
        return self.readonly_fields

    def save_model(self, request, obj, form, change): 
        obj.created_by= request.user
        obj.save()
        if obj.pinned_custom_field.all().count() == 0 and obj.schemaform:
            data = json.loads(form.cleaned_data["schemaform"])["form"]
            for position, field in enumerate(data):
                PinnedCustomField.objects.create(allowed_values=field["allowed_values"],
                                                custom_field_config=obj,
                                                field_type=field["field_type"],
                                                position=field["position"],
                                                name=field["key"],
                                                description=field["placeholder"])
        
                                                
                


    def log_change(self, request, object, message):
        """
        Log that an object has been successfully changed.
        The default implementation creates an admin LogEntry object.
        """
        super(CustomFieldConfigAdmin, self).log_change(request, object, message)
        cfr = ProjectResource()
        if object.__class__.__name__ == "CustomFieldConfig":
            schemaform = json.dumps(cfr.get_schema_form(object,"" ))
            object.schemaform = schemaform
            object.save()


    formfield_overrides = {
        models.CharField: {'widget': TextInput(attrs={'size':'20'})},
    }


class ProjectAdmin(ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'
    exclude= ["created_by"]
    actions = ['reindex']

    def save_model(self, request, obj, form, change): 
        obj.created_by = request.user
        obj.save()

    def reindex(self, request, queryset):
        cbr = CBHCompoundBatchResource()
        cbr.reindex_elasticsearch(request)
        self.message_user(request, "Successfully reindexed ChemReg compounds")
    reindex.short_description = "Reindex all compounds"



admin.site.register(CustomFieldConfig, CustomFieldConfigAdmin)
admin.site.register(Project, ProjectAdmin)
admin.site.register(SkinningConfig, SingletonModelAdmin)
