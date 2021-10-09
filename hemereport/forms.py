from django import forms
from .models import AllAbberations

class AllAbberationsForm(forms.ModelForm):
    class Meta:
        model = AllAbberations
        fields = ['genes','variants','tier','variant_type','vaf','coverage','transcript','locus','exon','length_bp']

        widgets = {
            'genes': forms.TextInput(attrs={ 'class': 'form-control'}),
            'variants': forms.TextInput(attrs={ 'class': 'form-control'}),
            'tier': forms.TextInput(attrs={ 'class': 'form-control'}),
            'variant_type': forms.TextInput(attrs={ 'class': 'form-control'}),
            'vaf': forms.TextInput(attrs={ 'class': 'form-control'}),
            'coverage': forms.TextInput(attrs={ 'class': 'form-control'}),
            'transcript': forms.TextInput(attrs={ 'class': 'form-control'}),
            'locus': forms.TextInput(attrs={ 'class': 'form-control'}),
            'exon': forms.TextInput(attrs={ 'class': 'form-control'}),
            'length_bp': forms.TextInput(attrs={ 'class': 'form-control'}),
            'variant_pmkb': forms.TextInput(attrs={ 'class': 'form-control'}),
        } 