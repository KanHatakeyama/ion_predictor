from django.contrib import admin
from import_export.admin import ImportExportModelAdmin
from . import models
from admin_numeric_filter.admin import RangeNumericFilter
from django_pandas.io import read_frame
from django.http import HttpResponse
from .ion_predictor.django_wrapper import auto_predictor

# Register your models here.
# class ChemicalAdmin(NestedModelAdmin):
class ChemicalAdmin(ImportExportModelAdmin):

    list_display = ["title", 
                    #"unique_name",
                      "subtitle", "smiles_thumbnail",
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


class CompositeAdmin(ImportExportModelAdmin):

    list_display = ["title",
                    # "unique_name",
                    #  "subtitle", 
                        "temperature",
                        "conductivity",
                        #"log_sigma",

                        "component1",
                        "component2",
                        "component3",
                        "mol_ratio",
                        "wt_ratio",
                    "created_at", "updated_at",
                    "tag_names",
                   ]
    list_filter = ["tags__name",  
                        ("temperature",RangeNumericFilter),
                        ("conductivity",RangeNumericFilter),
                   "created_at", "updated_at"]
    search_fields = ['title', "subtitle"]

    ordering = ["title"]
    save_as = True

    # show tags
    def tag_names(self, obj):
        return "\n".join([p.name for p in obj.tags.all()])

    actions = ['predict_conductivity']   

    def predict_conductivity(self, request, queryset):  

        composite_df=read_frame(queryset.all())
        compound_df=read_frame(models.Chemical.objects.all())
        predicted_df=auto_predictor.predict(composite_df,compound_df)

        #return csv
        response = HttpResponse(content_type='text/csv; charset=utf8')
        response['Content-Disposition'] = 'attachment; filename=users.csv'
        predicted_df.to_csv(path_or_buf=response, encoding='utf_8_sig', index=None)
        return response 

    predict_conductivity.short_description = 'predict conductivity'  

#add
admin.site.register(models.Chemical, ChemicalAdmin)
admin.site.register(models.Composite, CompositeAdmin)
admin.site.register(models.Tag)