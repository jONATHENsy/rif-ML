import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sentence_transformers import SentenceTransformer

# === 1. Load your CSV file (exported from ASReview) ===
# Replace this with your own CSV path
file_path = input("path to ur csv file ")
df = pd.read_csv(file_path, encoding='ISO-8859-1')

# === 2. Combine title and abstract into one text field ===
df['text'] = df[['Title', 'Abstract']].fillna('').agg(' '.join, axis=1)

# === 3. Extract labeled data (included == 1 or 0) as training set ===
df_labeled = df[df['included'].isin([0, 1])]
X_train_text = df_labeled['text'].tolist()
y_train = df_labeled['included']

# === 4. Extract unlabeled data (included is NaN) as test set ===
df_unlabeled = df[df['included'].isna()].copy()
X_unlabeled_text = df_unlabeled['text'].tolist()

# === 5. Generate sentence embeddings using Sentence-BERT ===
print("start sentence embeddings.....")
embedder = SentenceTransformer('all-MiniLM-L6-v2')  # You can try other models like 'all-mpnet-base-v2'
X_train_embed = embedder.encode(X_train_text, show_progress_bar=True)
X_unlabeled_embed = embedder.encode(X_unlabeled_text, show_progress_bar=True)

# === 6. Normalize embeddings for better classifier performance ===
scaler = StandardScaler()
X_train_embed = scaler.fit_transform(X_train_embed)
X_unlabeled_embed = scaler.transform(X_unlabeled_embed)

# === 7. Train Logistic Regression model ===
clf = LogisticRegression(max_iter=1000, class_weight='balanced')
clf.fit(X_train_embed, y_train)

# === 8. Predict relevance of unlabeled papers ===
probs = clf.predict_proba(X_unlabeled_embed)[:, 1]
df_unlabeled['predicted_relevance'] = probs
df_unlabeled['predicted_label'] = (probs >= 0.58).astype(int)  # You can adjust this threshold

# === 9. Sort and export results ===
df_result = df_unlabeled.sort_values(by='predicted_relevance', ascending=False)

# Choose which columns to export
columns_to_export = [
    'record_id', 'Authors', 'Title', 'Source Title', 'Abstract',
    'Publication Year', 'DOI',
    'UT (Unique WOS ID)', 'predicted_relevance', 'predicted_label'
]

# === Output to CSV ===
output_path = input("Enter the path to save the prediction result CSV (e.g., predicted_relevant.csv): ")
df_result[columns_to_export].to_csv(output_path, index=False)

print(f"✅ Done! saved to: {output_path}")
