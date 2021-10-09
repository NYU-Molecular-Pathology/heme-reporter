from django.contrib import admin
from .models import PMKBdb, AllAbberations, VariantTiering, VariantSpecificComment

# Register your models here.
admin.site.register(PMKBdb)
admin.site.register(AllAbberations)
admin.site.register(VariantTiering)
admin.site.register(VariantSpecificComment)

