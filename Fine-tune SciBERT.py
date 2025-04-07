import pandas as pd
import numpy as np
import torch
from datasets import Dataset
from transformers import AutoTokenizer, AutoModelForSequenceClassification, TrainingArguments, Trainer
from sklearn.metrics import accuracy_score, precision_recall_fscore_support

# 1. Load ASReview-exported dataset
file_path = r"D:\2025s1\BIOX7011\dlselect\asreview_dataset_all_Recovered Project.csv"
df = pd.read_csv(file_path, encoding='ISO-8859-1')

# 2. Combine Title + Abstract into a text field
df['text'] = df[['Title', 'Abstract']].fillna('').agg(' '.join, axis=1)

# 3. Prepare labeled training data
df_labeled = df[df['included'].isin([0, 1])].copy()
df_labeled = df_labeled[['text', 'included']].rename(columns={'included': 'label'})

df_labeled['label'] = df_labeled['label'].astype(int)

# 4. Prepare unlabeled data for prediction later
df_unlabeled = df[df['included'].isna()].copy()
df_unlabeled['text'] = df_unlabeled[['Title', 'Abstract']].fillna('').agg(' '.join, axis=1)

# 5. Load SciBERT tokenizer and model
model_name = "allenai/scibert_scivocab_uncased"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)

# 6. Tokenize the data
def preprocess_function(examples):
    return tokenizer(examples["text"], truncation=True, padding="max_length", max_length=256)

dataset = Dataset.from_pandas(df_labeled)
tokenized_dataset = dataset.map(preprocess_function, batched=True)

# 7. Train/test split (you can also use full set to train)
split_dataset = tokenized_dataset.train_test_split(test_size=0.1)

# 8. Metrics function
def compute_metrics(eval_pred):
    logits, labels = eval_pred
    preds = np.argmax(logits, axis=1)
    precision, recall, f1, _ = precision_recall_fscore_support(labels, preds, average="binary")
    acc = accuracy_score(labels, preds)
    return {"accuracy": acc, "precision": precision, "recall": recall, "f1": f1}

# 9. Training arguments
training_args = TrainingArguments(
    output_dir="./scibert_output",
    learning_rate=2e-5,
    per_device_train_batch_size=8,
    per_device_eval_batch_size=8,
    num_train_epochs=4,
    evaluation_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
    weight_decay=0.01,
    logging_dir="./logs",
)

# 10. Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=split_dataset["train"],
    eval_dataset=split_dataset["test"],
    tokenizer=tokenizer,
    compute_metrics=compute_metrics,
)

# 11. Train!
trainer.train()

# 12. Predict relevance on unlabeled data
tokenized_unlabeled = Dataset.from_pandas(df_unlabeled[['text']])
tokenized_unlabeled = tokenized_unlabeled.map(preprocess_function, batched=True)

predictions = trainer.predict(tokenized_unlabeled)
probs = torch.nn.functional.softmax(torch.tensor(predictions.predictions), dim=1)[:, 1].numpy()

df_unlabeled['predicted_relevance'] = probs
df_unlabeled['predicted_label'] = (df_unlabeled['predicted_relevance'] >= 0.52).astype(int)

# 13. Sort and export
df_result = df_unlabeled.sort_values(by='predicted_relevance', ascending=False)

columns_to_export = [
    'record_id', 'Authors', 'Title', 'Source Title', 'Abstract',
    'Publication Year', 'DOI',
    'UT (Unique WOS ID)', 'predicted_relevance', 'predicted_label'
]

output_path = r"D:\2025s1\BIOX7011\dlselect\aftermannual_finetuned_scibert_predictions.csv"
df_result[columns_to_export].to_csv(output_path, index=False)

print(f"✅ Done! Fine-tuned SciBERT predictions saved to:\n📁 {output_path}")
