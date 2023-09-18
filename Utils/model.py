import torch
import torch.nn as nn
# from transformers import BertModel, BertTokenizer
import math

class BiLSTM_Attention(nn.Module):
    def __init__(self, input_dim=64, hidden_dim=128, num_layers=3, output_dim=2, dropout_prob=0.5):
        super(BiLSTM_Attention, self).__init__()

        # Define fc
        self.linear = nn.Linear(27, input_dim)

        # Define the LSTM layer
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers, bidirectional=True, batch_first=True)

        # Define the attention layer
        self.attention = nn.Linear(hidden_dim * 2, 1)

        # Define the dropout layer
        self.dropout = nn.Dropout(dropout_prob)

        # Define the fully connected layer
        self.fc = nn.Linear(hidden_dim * 2, output_dim)

    def forward(self, x):
        # Get the length of the input sequence
        seq_len = x.shape[1]

        x = self.linear(x)
        # Pass the input through the LSTM layer
        lstm_out, _ = self.lstm(x)

        # Compute the attention weights
        attention_scores = self.attention(lstm_out).squeeze()
        attention_weights = nn.functional.softmax(attention_scores, dim=1).unsqueeze(2)

        # Compute the weighted sum of the output of the LSTM layer
        weighted_sum = torch.bmm(lstm_out.transpose(1, 2), attention_weights).squeeze()

        # Apply dropout
        out = self.dropout(weighted_sum)

        # Pass the output through a fully connected layer
        out = self.fc(out)

        return out


# Define transformer model
class TransformerClassifier(nn.Module):
    def __init__(self, feature_dim=27, num_classes=2, max_seq_len=6, num_layers=1, hidden_dim=128, num_heads=8, dropout=0.5):
        super(TransformerClassifier, self).__init__()
        self.embedding = nn.Linear(feature_dim, hidden_dim)
        self.positional_encodings = self.generate_positional_encodings(max_seq_len, hidden_dim).to(device)
        self.transformer_layers = nn.TransformerEncoderLayer(hidden_dim, nhead=num_heads)
        self.transformer = nn.TransformerEncoder(self.transformer_layers, num_layers=num_layers)
        self.fc = nn.Linear(hidden_dim * max_seq_len, num_classes)

    def generate_positional_encodings(self, max_seq_len, hidden_dim):
        position = torch.arange(0, max_seq_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, hidden_dim, 2) * -(math.log(10000.0) / hidden_dim))
        positional_encodings = torch.zeros((max_seq_len, hidden_dim))
        positional_encodings[:, 0::2] = torch.sin(position * div_term)
        positional_encodings[:, 1::2] = torch.cos(position * div_term)
        return positional_encodings.unsqueeze(0)

    def forward(self, x):
        batch_size, seq_len, feature_dim = x.size()
        x = self.embedding(x) + self.positional_encodings[:, :seq_len, :]
        x = x.permute(1, 0, 2)
        output = self.transformer(x)
        output = output.permute(1, 0, 2).reshape(batch_size, -1)
        logits = self.fc(output)
        return logits

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")

def load_checkpoint(load_path, model, optimizer):
    if load_path == None:
        return

    state_dict = torch.load(load_path, map_location=device)
    print(f'Model loaded from <== {load_path}')

    model.load_state_dict(state_dict['model_state_dict'])
    optimizer.load_state_dict(state_dict['optimizer_state_dict'])

    return state_dict['valid_loss']
