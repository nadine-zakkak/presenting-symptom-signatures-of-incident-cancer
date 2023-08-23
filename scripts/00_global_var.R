# Nadine Zakkak
# Script that contains global variables required for analysis and output of figures and tables

# Groups Ordering ----
# order based on heat map visual
cancergrp_order <- c("Head and neck", "Upper GI", "Lower GI", "HPB", "Respiratory",
                     "Urological", "Haematological", "CNS", "Sarcoma", "Skin",
                     "Ocular", "Breast", "Gynaecological", "Prostate and other male organs", 
                     "Other malignant neoplasms", "Unknown primary")
cancer_order <- c("Larynx", "Oral cavity", "Oropharynx", "Thyroid", "Other head and neck", 
                  "Oesophagus", "Stomach", 
                  "Anal", "Colon", "Rectum", "Small intestine", 
                  "Liver", "Pancreas", "Other HPB",
                  "Lung", "Mesothelioma", 
                  "Bladder", "Kidney", "Ureteric and other urinary", 
                  "Acute leukaemia", "Chronic lymphocytic leukaemia", "Hodgkin lymphoma", "Multiple myeloma", "Non-Hodgkin lymphoma", "Other haematological", 
                  "CNS", 
                  "Bone sarcoma", "Connective and soft tissue sarcoma", 
                  "Melanoma", 
                  "Ocular", 
                  "Breast", 
                  "Cervix", "Ovary", "Uterus", "Vulva/Vagina", 
                  "Penile", "Prostate", "Testicular",
                  "Other malignant neoplasms","Unknown primary") #enforcing "Other XXX" to be at end of each group

sxgrp_order <- c("Non-specific", "Lump/mass/lymph node", "Ulceration",
                 "Upper abdominal", "Lower abdominal",
                 "Respiratory", "Urological", "Central nervous system", "Musculoskeletal",
                 "Skin Lesion", "Breast Symptoms", "Female specific", "Male specific",
                 "None recorded", "summary")
