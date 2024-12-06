#!/usr/bin/env python

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import lightning as L

import os, argparse

# Added this class becase datautils no longer exists
class MyTrainDataset(Dataset):
    def __init__(self, size):
        self.size = size
        self.x = torch.rand(size,20)
        self.y = torch.rand(size,1)
    def __len__(self):
        return self.size
    def __getitem__(self, index):
        return (self.x[index,:], self.y[index])

class LanguageModel(L.LightningModule):
    def __init__(self):
        super().__init__()
        layers = [torch.nn.Linear(20, 512),torch.nn.Linear(512,2048)]
        layers += [torch.nn.Linear(2048,2048) for i in range(13)]
        layers.append(torch.nn.Linear(2048, 1))
        self.model = torch.nn.Sequential(*layers)

    def training_step(self, batch, batch_idx):
        input, target = batch
        output = self.model(input)
        loss = F.cross_entropy(output, target)
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def configure_optimizers(self):
        return torch.optim.SGD(self.model.parameters(), lr=1e-3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='simple distributed training job')
    parser.add_argument('--epochs', default=3, type=int, help='Total epochs to train the model')
    parser.add_argument('--batch_size', default=512, type=int, help='Input batch size on each device (default: 512)')
    parser.add_argument('--data_size', default=2048*64, type=int, help='Number of items in dataset (default: %(default)s)')
    parser.add_argument('-N', default=1, type=int, help='Number of nodes')
    parser.add_argument('-p', default=1, type=int, help='Number of gpus per node')

    args = parser.parse_args()
    
    train_dataloader = DataLoader(MyTrainDataset(args.data_size), batch_size=args.batch_size)
    model = LanguageModel()

    # Trainer
    trainer = L.Trainer(max_epochs=args.epochs, devices=args.p, accelerator="gpu", num_nodes=args.N)
    trainer.fit(model, train_dataloader)
