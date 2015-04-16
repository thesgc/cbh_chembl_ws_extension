from django.contrib import admin
from cbh_chembl_model_extension.models import Project, PinnedCustomField, CustomFieldConfig

from django.contrib.admin import ModelAdmin


from django.forms.widgets import HiddenInput, TextInput
from django.db import models
import json

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


#Make a template have to be chosen in order to create a schema and make it impossible to edit schemas once created then versioning not needed



class CustomFieldConfigAdmin(ModelAdmin):

    exclude= ["created_by", ]

    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created' 
    inlines = [PinnedCustomFieldInline,]
    def save_model(self, request, obj, form, change): 
        obj.created_by = request.user
        print "test"
        obj.save()
    formfield_overrides = {
        models.CharField: {'widget': TextInput(attrs={'size':'20'})},
    }



    # def save_related(self, request, form, formsets, change):
    #     """
    #     Given the ``HttpRequest``, the parent ``ModelForm`` instance, the
    #     list of inline formsets and a boolean value based on whether the
    #     parent is being added or changed, save the related objects to the
    #     database. Note that at this point save_form() and save_model() have
    #     already been called.
    #     """
    #     form.save_m2m()
    #     for formset in formsets:
    #         instances = self.save_formset(request, form, formset, change=change)
    #         fields = []
    #         for f in formset.queryset:
    #             f.field_key = slugify(f.name).replace("-", "_")
    #             f.save()
    #             fields.append(f.get_fields())
    #         schemaform = {
    #                         "schema" :{
    #                                     "type" : "object",
    #                                     "properties"   :  dict((field[0],field[1]) for field in fields),
    #                                     "required" : [field[0] for field in fields if field[2] is True]
    #                         },
    #                         "form" : [field[0] for field in fields]
    #                     }
    #         formset.instance.schemaform = json.dumps(schemaform)
    #         formset.instance.save()
    #         break # There is only one formset




     # get_template('templates/email.html').render(
     #        Context()
     #    ),

     
        # javascript = """var schema = 
        # if data["format"] == UISELECT:
        #     data["options"] =  {
        #           "tagging": "tagFunction" ,
        #           "taggingLabel": "(adding new)",
        #           "taggingTokens": "",
        #        },




class ProjectAdmin(ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'
    exclude= ["created_by"]

    def save_model(self, request, obj, form, change): 
        obj.created_by = request.user
        obj.save()


admin.site.register(CustomFieldConfig, CustomFieldConfigAdmin)
admin.site.register(Project, ProjectAdmin)
