from django.contrib import admin
from import_export.admin import ImportExportModelAdmin
from . import models
from django.forms import Textarea

# Register your models here.
# class ChemicalAdmin(NestedModelAdmin):
class ChemicalAdmin(ImportExportModelAdmin):

    list_display = ["title", "unique_name", "subtitle", "smiles_thumbnail",
                    "created_at", "updated_at",
                    "tag_names",
                   ]
    list_filter = ["tags__name",  
                   "created_at", "updated_at"]
    search_fields = ['title', "subtitle", "SMILES"]

    readonly_fields = ('smiles_preview',"smiles_info",)
    ordering = ["title"]
    save_as = True

    #actions = ["show_table", 'show_json']

    #formfield_overrides = {
    #        models.TextField: {'widget': Textarea(attrs={'rows':1, 'cols':100})},
    #    }


    # show tags
    def tag_names(self, obj):
        return "\n".join([p.name for p in obj.tags.all()])

    # show chemical structure image
    def smiles_preview(self, obj):
        return obj.smiles_preview
    smiles_preview.short_description = 'Structure'
    smiles_preview.allow_tags = True

    def smiles_thumbnail(self, obj):
        return obj.smiles_thumbnail
    smiles_thumbnail.short_description = 'Structure'
    smiles_thumbnail.allow_tags = True

    # show unique name
    def unique_name(self, obj):
        return obj.unique_name
    unique_name.short_description = 'unique_name'
    unique_name.allow_tags = True

    def smiles_info(self, obj):
        return obj.smiles_info
    smiles_info.short_description = 'smiles_info'
    smiles_info.allow_tags = True


#add
admin.site.register(models.Chemical, ChemicalAdmin)
admin.site.register(models.Tag)