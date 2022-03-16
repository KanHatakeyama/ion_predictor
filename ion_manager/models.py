from django.db import models
#from ckeditor_uploader.fields import RichTextUploadingField
from .ion_predictor.django_wrapper import parse_smiles
from django.utils.html import mark_safe
import numpy as np

# Create your models here.
# tag
class Tag(models.Model):
    name = models.CharField(max_length=32)

    def __str__(self):
        return self.name



# chemical
class Chemical(models.Model):
    title = models.CharField(max_length=200,unique=True)
    subtitle = models.CharField(max_length=200, null=True, blank=True)
    tags = models.ManyToManyField(Tag, null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    obtained_date = models.DateField(null=True, blank=True)
    #made_by = models.CharField(max_length=400, null=True, blank=True)
    SMILES = models.TextField(max_length=4000, null=True, blank=True)
    mol_ratio= models.CharField(max_length=400, null=True, blank=True)
    wt_ratio= models.CharField(max_length=400, null=True, blank=True)
    Mw= models.FloatField(max_length=400, null=True, blank=True)
    Mn= models.FloatField(max_length=400, null=True, blank=True)
    MwMn= models.FloatField(max_length=400, null=True, blank=True)
    Polymn_deg=models.FloatField(max_length=400, null=True, blank=True)
    Structure= models.CharField(max_length=400, null=True, blank=True)
    melting_temp=models.CharField(max_length=400, null=True, blank=True)
    Tg=models.CharField(max_length=400, null=True, blank=True)

    #company = models.CharField(max_length=200, null=True, blank=True)
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



# composite 
class Composite(models.Model):
    title = models.CharField(max_length=200,unique=True)
    subtitle = models.CharField(max_length=200, null=True, blank=True)
    tags = models.ManyToManyField(Tag, null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    obtained_date = models.DateField(null=True, blank=True)


    temperature=models.FloatField(max_length=400)
    conductivity=models.FloatField(max_length=400)

    component1 = models.ForeignKey(Chemical,models.SET_NULL,blank=True,null=True,related_name="chem_component1")
    component2 = models.ForeignKey(Chemical,models.SET_NULL,blank=True,null=True,related_name="chem_component2")
    component3 = models.ForeignKey(Chemical,models.SET_NULL,blank=True,null=True,related_name="chem_component3")
    component4 = models.ForeignKey(Chemical,models.SET_NULL,blank=True,null=True,related_name="chem_component4")
    component5 = models.ForeignKey(Chemical,models.SET_NULL,blank=True,null=True,related_name="chem_component5")
    component6 = models.ForeignKey(Chemical,models.SET_NULL,blank=True,null=True,related_name="chem_component6")


    mol_ratio= models.CharField(max_length=400, null=True, blank=True)
    wt_ratio= models.CharField(max_length=400, null=True, blank=True)

    melting_temp=models.CharField(max_length=400, null=True, blank=True)
    inorg_name=models.CharField(max_length=400, null=True, blank=True)
    inorg_contain_ratio=models.FloatField(max_length=400, null=True, blank=True)
    crystallinity=models.FloatField(max_length=400, null=True, blank=True)
    Tg=models.FloatField(max_length=400, null=True, blank=True)
    mp=models.FloatField(max_length=400, null=True, blank=True)



    reference = models.CharField(max_length=400, null=True, blank=True)
    special_memo = models.TextField(max_length=4000,blank=True, null=True)

    class Meta:
        verbose_name = "Composite"

    def __str__(self) -> str:
        return self.title


    @property
    def unique_name(self):
        return str(self.pk)+"_"+str(self.title)

    @property
    def temperature_(self):
        return round(self.temperature,1)

    @property
    def log_sigma(self):
        return round(np.log10(self.conductivity),2)