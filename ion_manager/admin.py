from django.contrib import admin
from import_export.admin import ImportExportModelAdmin
from . import models
from admin_numeric_filter.admin import RangeNumericFilter
from django_pandas.io import read_frame
from django.http import HttpResponse
from .ion_predictor.django_wrapper import auto_predictor


def df_to_csv(df, filename="users.csv"):
    response = HttpResponse(content_type='text/csv; charset=utf8')
    response['Content-Disposition'] = f'attachment; filename={filename}'
    df.to_csv(path_or_buf=response, encoding='utf_8_sig', index=None)
    return response


# Register your models here.
# class ChemicalAdmin(NestedModelAdmin):
class ChemicalAdmin(ImportExportModelAdmin):

    list_display = ["title",
                    "tag_names",
                    "subtitle", "smiles_thumbnail",
                    "created_at", "updated_at",
                    ]
    list_filter = ["tags__name",
                   "created_at", "updated_at"]
    search_fields = ['title', "subtitle", "SMILES"]

    readonly_fields = ('smiles_preview', "smiles_info",)
    ordering = ["-updated_at", "pk"]
    save_as = True

    actions = ["dump"]

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

    def dump(self, request, queryset):
        df = read_frame(queryset.all())
        df = auto_predictor.pre_convert_compound(df)
        return df_to_csv(df, filename="compounds.csv")
    dump.short_description = 'dump for ML'


class CompositeAdmin(ImportExportModelAdmin):

    list_display = ["title",
                    "tag_names",
                    "temperature",
                    "conductivity",
                    "component1",
                    "component2",
                    "component3",
                    "mol_ratio",
                    "wt_ratio",
                    "created_at", "updated_at",
                    ]
    list_filter = ["tags__name",
                   ("temperature", RangeNumericFilter),
                   ("conductivity", RangeNumericFilter),
                   "created_at", "updated_at"]
    search_fields = ['title', "subtitle"]

    ordering = ["-updated_at", "pk"]
    save_as = True

    # show tags
    def tag_names(self, obj):
        return "\n".join([p.name for p in obj.tags.all()])

    actions = ['predict_conductivity', "dump"]

    def predict_conductivity(self, request, queryset):

        composite_df = read_frame(queryset.all())
        compound_df = read_frame(models.Chemical.objects.all())
        predicted_df = auto_predictor.predict(composite_df, compound_df)

        # return csv
        return df_to_csv(predicted_df, filename="prediction.csv")

    predict_conductivity.short_description = 'predict conductivity'

    def dump(self, request, queryset):
        composite_df = read_frame(queryset.all())
        composite_df = auto_predictor.pre_convert_composite(composite_df)
        return df_to_csv(composite_df, filename="composite.csv")
    dump.short_description = 'dump for ML'


# add
admin.site.register(models.Chemical, ChemicalAdmin)
admin.site.register(models.Composite, CompositeAdmin)
admin.site.register(models.Tag)
