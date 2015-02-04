from django.contrib import admin
from cbh_chembl_model_extension.models import Project, PinnedCustomField

from django.contrib.admin import ModelAdmin

class PinnedCustomFieldInline(admin.TabularInline):
    model = PinnedCustomField


class ProjectAdmin(ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'
    inlines = [PinnedCustomFieldInline,]

admin.site.register(Project, ProjectAdmin)
