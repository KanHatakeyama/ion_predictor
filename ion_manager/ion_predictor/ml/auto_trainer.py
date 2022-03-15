
from rdkit import Chem
import joblib
from sklearn.model_selection import train_test_split
import torch.nn as nn
import torch.optim as optim
import torch
from torch.utils.data import DataLoader
import dgl
from tqdm import tqdm
import matplotlib.pyplot as plt
import yaml

from .data_converter import mol2dgl_single
from .normalGNN import collate,Regressor,ATOM_FDIM


def load_trained_model(settings):

    #load model
    model = Regressor(ATOM_FDIM,256,32, int(settings["pretrain_y_dim"]))
    model.load_state_dict(torch.load(settings["model_path"], map_location=torch.device('cpu')))
    model.eval()
    model.hidden_mode=True

    return model

def auto_prepare_model(settings):

    model_path=settings["model_path"]

    train_graphs, test_graphs, train_y, test_y=prepare_dataset(settings)
    model = Regressor(ATOM_FDIM,256,32, settings["pretrain_y_dim"])

    model,epoch_losses=train_gnn_model(train_graphs,
                        train_y,
                        model,
                        settings
    )

    #save
    torch.save(model.to('cpu').state_dict(), model_path)
    plt.plot(epoch_losses, c='b')

def prepare_dataset(settings):
    """
    prepare dataset for ML

    Parameters
    ----------
    smiles_path: str
        path of SMILES data used for pretraining (str list)
        
    descriptor_path: str
        path of descriptor data (np array)
    
    test_size: float
        test size (use train_test_split func of sklearn)

    random_state: int
        random state

    Returns
    -------
    train_graphs, test_graphs, train_y, test_y: DGL graph, DGL graph, np array, np array
        training and testing data

    """

    descriptor_path=settings["descriptor_path"]
    smiles_path=settings["smiles_path"]
    test_size=settings["pretrain_test_size"]
    random_state=settings["random_state"]

    
    #X (DGL graph data)
    smiles_list=joblib.load(smiles_path)
    mols=[Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    graphs= mol2dgl_single(mols)

    #y (molecular descriptor)
    desc_array=joblib.load(descriptor_path)
    desc_array=desc_array.astype("float32")

    #split
    train_graphs, test_graphs, train_y, test_y= train_test_split(
    graphs, desc_array, test_size=test_size, random_state=random_state)

    return train_graphs, test_graphs, train_y, test_y



def train_gnn_model(train_graphs,
                     train_y,
                     model,
                     settings,
                     loss_func = nn.MSELoss(),
           ):

    model_path =settings["model_path"]

    #prepare dataloader
    dataset = list(zip(train_graphs, train_y))
    data_loader = DataLoader(dataset, batch_size=settings["batch_size"], shuffle=True, collate_fn=collate)

    #prepare regressor
    optimizer = optim.Adam(model.parameters(), lr=settings["lr"])
    model.train()

    #train epochs
    epoch_losses = []
    for epoch in tqdm(range(settings["epochs"])):
        epoch_loss = 0
        for i, (bg, label) in enumerate(data_loader):
            bg.set_e_initializer(dgl.init.zero_initializer)
            bg.set_n_initializer(dgl.init.zero_initializer)        
            pred = model(bg)
            loss = loss_func(pred, label)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.detach().item()
        epoch_loss /= (i + 1)
        if (epoch+1) % 20 == 0:
            print('Epoch {}, loss {:.4f}'.format(epoch+1, epoch_loss))
        epoch_losses.append(epoch_loss)

    return model,epoch_losses


"""
# test data check 
model.eval()
test_bg = dgl.batch(test_graphs)
test_y_tensor = torch.tensor(test_y).float().view(-1,1)
test_bg.set_e_initializer(dgl.init.zero_initializer)
test_bg.set_n_initializer(dgl.init.zero_initializer)
pred_y = model(test_bg).detach().numpy()
plt.plot(test_y[:,0],pred_y[:,0],"o")

"""