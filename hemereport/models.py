from django.db import models
from django.conf import settings

# Create your models here.
class PMKBdb(models.Model):
    gene = models.CharField(blank=False, max_length=255)
    tumor_type = models.TextField(blank=False)
    tissue_type = models.CharField(blank=False, max_length=255)
    variant = models.CharField(blank=False, max_length=255)
    tier = models.IntegerField()
    interpretations = models.TextField(blank=False)
    citations = models.TextField(blank=False)
    
    def __str__(self):
        return self.gene

class AllAbberations(models.Model):
    runid = models.CharField(blank=False,default='',max_length=255)
    sample = models.CharField(blank=False, max_length=255)
    genes = models.CharField(blank=False, max_length=255)
    variant_type = models.CharField(blank=False, default='NA',max_length=255)
    transcript = models.CharField(default='NA',max_length=255)
    variants = models.CharField(default='NA',max_length=255)
    length_bp = models.IntegerField(default=0)
    vaf = models.DecimalField(max_digits=6,decimal_places=2,default=0)
    exon = models.IntegerField(default=0)
    amino_acid_change = models.CharField(default='NA',max_length=255)
    coverage = models.IntegerField(default=0)
    tier = models.IntegerField()
    locus = models.CharField(blank=False, max_length=255)
    interpretations = models.TextField(default='NA')
    citations = models.TextField(default='NA')
    variant_pmkb = models.TextField(default='NA')
    class Meta:
        db_table="allabberations"
        
class VariantTiering(models.Model):
    abberation_id = models.IntegerField(blank=False,default='0')
    user = models.CharField(max_length=20)
    runid = models.CharField(blank=False,default='',max_length=255)
    sample = models.CharField(blank=False, max_length=255)
    genes = models.CharField(blank=False, max_length=255)
    variants = models.CharField(default='NA',max_length=255)
    tier = models.IntegerField()
    
class VariantSpecificComment(models.Model):
    genes = models.CharField(blank=False, max_length=255)
    variants = models.CharField(default='NA',max_length=255)
    tier = models.IntegerField()
    interpretations = models.TextField(default='NA')
    citations = models.TextField(default='NA')