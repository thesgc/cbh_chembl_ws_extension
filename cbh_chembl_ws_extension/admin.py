# -*- coding: utf-8 -*-

from django.contrib import admin
from cbh_core_model.models import Project, PinnedCustomField, CustomFieldConfig, SkinningConfig, ProjectType

from django.contrib.admin import ModelAdmin
from cbh_chembl_ws_extension.projects import ChemregProjectResource
from cbh_chembl_ws_extension.compounds import CBHCompoundBatchResource

from django.forms.widgets import HiddenInput, TextInput
from django.db import models
import json
from solo.admin import SingletonModelAdmin

class ChemregProject(Project):
    class Meta:
        proxy = True


class ProjectAdmin(ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created', 'project_type')
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




admin.site.register(ChemregProject, ProjectAdmin)

