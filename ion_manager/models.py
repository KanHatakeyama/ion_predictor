from django.db import models
#from ckeditor_uploader.fields import RichTextUploadingField
from .ion_predictor.django_wrapper import parse_smiles
from django.utils.html import mark_safe

# Create your models here.
# tag
class Tag(models.Model):
    name = models.CharField(max_length=32)

    def __str__(self):
        return self.name


# chemical
class Chemical(models.Model):
    title = models.CharField(max_length=200)
    subtitle = models.CharField(max_length=200, null=True, blank=True)
    tags = models.ManyToManyField(Tag, null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    obtained_date = models.DateField(null=True, blank=True)
    made_by = models.CharField(max_length=400, null=True, blank=True)
    SMILES = models.TextField(max_length=4000, null=True, blank=True)
    company = models.CharField(max_length=200, null=True, blank=True)
    reference = models.CharField(max_length=400, null=True, blank=True)
    #special_memo = RichTextUploadingField(blank=True, null=True)
    special_memo = models.TextField(max_length=4000,blank=True, null=True)

    class Meta:
        verbose_name = "Chemical"

    def __str__(self) -> str:
        return self.title

        # this somehow causes error during integrating graphs (utils/experiment_utilities)
        # return str(self.pk)+"_"+str(self.title)

    # show chemical structures
    def prepare_image(self, smiles, size=None):
        try:
            if size is None:
                img = parse_smiles.smiles_to_buffer_img(smiles)
            else:
                img = parse_smiles.smiles_to_buffer_img(smiles, size=size)
            tag = "<img src='data:image/png;base64,{}'/>".format(img)
        except:
            tag = "<p>error parsing smiles</p>"
        return tag

    # polymer data parsing
    @property
    def smiles_info(self):
        return parse_smiles.smiles_info(self.SMILES)

    @property
    def smiles_preview(self):
        if self.SMILES:
            return mark_safe(self.prepare_image(self.SMILES))
        return ""

    @property
    def smiles_thumbnail(self):
        if self.SMILES:
            return mark_safe(self.prepare_image(self.SMILES, size=100))
        return ""

    @property
    def unique_name(self):
        return str(self.pk)+"_"+str(self.title)