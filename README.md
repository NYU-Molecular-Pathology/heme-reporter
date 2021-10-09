# heme-reporter

Clinical reporting application for Oncomine Myeloid assay results.


<img width="1634" alt="heme-reporter-home" src="https://user-images.githubusercontent.com/25512313/136666925-eaea424a-f484-4af2-8406-4d5a697eae32.png">

## Overview

The Thermo Fisher S5 Next-Gen Sequencing platform can be used for the identification of structural variants, mutations, and other abnormalities in DNA and RNA samples using on-board Ion Reporter (IR) data analysis software suite. However, Ion Reporter does not include clinical interpretation of variants, which is a desired feature for users of the platform.

Clinicians at NYU have compiled their own sets of interpretations and tiers of clinically relevant variants, for use with sign-out of cases using the S5 system.

**heme-reporter** app combines a pre-saved version of PMKB database along with custom NYU variant interpretations in order to match clinical interpretations with user-provided variants exported from the Ion-Reporter in .csv format. This web application allows users to upload their .csv formatted variant export file and find interpretations for each variant, matched on gene comparing to pre-saved PMKB data and helps them to edit, update, tier and preview variants in their desired reporting template which can then be exported to word document. This is then copied over to EPIC system for signing out clinical cases hence reducing the amount of time it takes to manually report clinical cases.
