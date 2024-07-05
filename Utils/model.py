import torch
import torch.nn as nn

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


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")

def load_checkpoint(load_path, model, optimizer):
    if load_path == None:
        return

    state_dict = torch.load(load_path, map_location=device)
    print(f'Model loaded from <== {load_path}')

    model.load_state_dict(state_dict['model_state_dict'])
    optimizer.load_state_dict(state_dict['optimizer_state_dict'])

    return state_dict['valid_loss']
