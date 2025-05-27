# scRNA-seq Streamlit App

Αυτή η εφαρμογή υλοποιεί μια διαδραστική πλατφόρμα ανάλυσης δεδομένων scRNA-seq (single-cell RNA sequencing) με χρήση της βιβλιοθήκης Scanpy και Streamlit.

## 🔬 Λειτουργίες
- Προεπεξεργασία δεδομένων (QC filtering)
- Ανάλυση clustering (Leiden)
- UMAP προβολή και οπτικοποιήσεις
- Violin plots για γονίδια ενδιαφέροντος
- Dockerized εκτέλεση για ευκολία εγκατάστασης

## 🐳 Εκτέλεση με Docker

```bash
docker-compose up
```

Ή χωρίς Docker (τοπικά):

```bash
streamlit run streamlit_pipeline_app.py
```

## 📁 Δεδομένα

⚠️ Το αρχείο `pancreas_data.h5ad` δεν περιλαμβάνεται στο GitHub λόγω μεγέθους (>100MB).  
Μπορείτε να το κατεβάσετε από εδώ:  
🔗 [Λήψη από Google Drive](https://drive.google.com/file/d/1lec5twDcoR_JNMrkId7czwsCuHLXHwax/view?usp=drive_link)

## 👥 Μέλη Ομάδας
Δείτε και το αρχείο `AUTHORS.txt` για αναλυτική συνεισφορά.
