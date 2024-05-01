#!/usr/bin/env python

# Adapted from https://github.com/pytorch/examples/blob/main/distributed/ddp-tutorial-series/multinode.py

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import torch.multiprocessing as mp
from torch.utils.data.distributed import DistributedSampler
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.distributed import init_process_group, destroy_process_group
import os

def ddp_setup():
    init_process_group(backend="nccl")
    torch.cuda.set_device(int(os.environ["LOCAL_RANK"]))

# Added this class becase datautils no longer exists
class MyTrainDataset(Dataset):
    def __init__(self, size):
        self.size = size
        self.data = [(torch.rand(20), torch.rand(1)) for _ in range(size)]
    def __len__(self):
        return self.size
    def __getitem__(self, index):
        return self.data[index]

class Trainer:
    def __init__(
        self,
        model: torch.nn.Module,
        train_data: DataLoader,
        optimizer: torch.optim.Optimizer,
        snapshot_path: str,
    ) -> None:
        self.local_rank = int(os.environ["LOCAL_RANK"])
        self.global_rank = int(os.environ["RANK"])
        self.model = model.to(self.local_rank)
        self.train_data = train_data
        self.optimizer = optimizer
        self.epochs_run = 0
        self.snapshot_path = snapshot_path
        self.model = DDP(self.model, device_ids=[self.local_rank])

    def _load_snapshot(self, snapshot_path):
        loc = f"cuda:{self.local_rank}"
        snapshot = torch.load(snapshot_path, map_location=loc)
        self.model.load_state_dict(snapshot["MODEL_STATE"])
        self.epochs_run = snapshot["EPOCHS_RUN"]
        print(f"Resuming training from snapshot at Epoch {self.epochs_run}")

    def _run_batch(self, source, targets):
        self.optimizer.zero_grad()
        output = self.model(source)
        loss = F.cross_entropy(output, targets)
        loss.backward()
        self.optimizer.step()

    def _run_epoch(self, epoch):
        b_sz = len(next(iter(self.train_data))[0])
        print(f"[GPU{self.global_rank}] Epoch {epoch} | Batchsize: {b_sz} | Steps: {len(self.train_data)}")
        self.train_data.sampler.set_epoch(epoch)
        for source, targets in self.train_data:
            source = source.to(self.local_rank)
            targets = targets.to(self.local_rank)
            self._run_batch(source, targets)

    def train(self, max_epochs: int):
        for epoch in range(self.epochs_run, max_epochs):
            self._run_epoch(epoch)

def load_train_objs():
    train_set = MyTrainDataset(2048*2048)  # load your dataset
    # Make layers for model
    layers = [torch.nn.Linear(20, 512),torch.nn.Linear(512,2048)]
    layers += [torch.nn.Linear(2048,2048) for i in range(13)]
    layers.append(torch.nn.Linear(2048, 1))
    model = torch.nn.Sequential(*layers)
    optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)
    return train_set, model, optimizer

def prepare_dataloader(dataset: Dataset, batch_size: int):
    return DataLoader(
        dataset,
        batch_size=batch_size,
        pin_memory=True,
        shuffle=False,
        # Scales steps with number of workers
        sampler=DistributedSampler(dataset)
    )

def main(total_epochs: int, batch_size: int, snapshot_path: str = "snapshot.pt"):
    ddp_setup()
    dataset, model, optimizer = load_train_objs()
    train_data = prepare_dataloader(dataset, batch_size)
    trainer = Trainer(model, train_data, optimizer, snapshot_path)
    trainer.train(total_epochs)
    destroy_process_group()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='simple distributed training job')
    parser.add_argument('--epochs', defaut=5, type=int, help='Total epochs to train the model')
    parser.add_argument('--batch_size', default=512, type=int, help='Input batch size on each device (default: 512)')
    args = parser.parse_args()
    
    main(args.epochs, args.batch_size)
