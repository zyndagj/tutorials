#!/usr/bin/env python

# Adapted from https://github.com/pytorch/examples/blob/main/distributed/ddp-tutorial-series/multinode.py

import torch, argparse, os
from time import time
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

from torch.utils.data.distributed import DistributedSampler
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.distributed import init_process_group, destroy_process_group

def ddp_setup():
    torch.cuda.set_device(int(os.environ["LOCAL_RANK"]))
    init_process_group(backend="nccl")

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

class Trainer:
    def __init__(
        self,
        model: torch.nn.Module,
        train_data: DataLoader,
        optimizer: torch.optim.Optimizer
    ) -> None:
        self.local_rank = int(os.environ["LOCAL_RANK"])
        self.global_rank = int(os.environ["RANK"])
        self.model = model.to(self.local_rank)
        self.train_data = train_data
        self.optimizer = optimizer
        self.epochs_run = 0
        self.model = DDP(self.model, device_ids=[self.local_rank])

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
            s = time()
            self._run_epoch(epoch)
            etime = time()-s
            if not self.global_rank:
                print(f"Finished epoch {epoch} in {etime:.2f} seconds")

def load_train_objs(data_size: int):
    train_set = MyTrainDataset(data_size)  # load your dataset
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

def main(total_epochs: int, batch_size: int, data_size: int):
    ddp_setup()
    gr = int(os.environ["RANK"])
    if not gr:
        print("Finished DDP setup")
        print(f"Dataset has {data_size} items")
    dataset, model, optimizer = load_train_objs(data_size)
    train_data = prepare_dataloader(dataset, batch_size)
    trainer = Trainer(model, train_data, optimizer)
    trainer.train(total_epochs)
    destroy_process_group()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='simple distributed training job')
    parser.add_argument('--epochs', default=3, type=int, help='Total epochs to train the model')
    parser.add_argument('--batch_size', default=512, type=int, help='Input batch size on each device (default: 512)')
    parser.add_argument('--data_size', default=2048*64, type=int, help='Number of items in dataset (default: %(default)s)')
    args = parser.parse_args()
    
    main(args.epochs, args.batch_size, args.data_size)
