import torch
import dgl
import dgl.function as fn
import torch.nn as nn
import torch.nn.functional as F

from .data_converter import ATOM_FDIM

"""
Convert molecular structures to 32 dim vectors by a neural net

"""
# Note that during graph decoding they don't predict stereochemistry-related
# characteristics (i.e. Chiral Atoms, E-Z, Cis-Trans).  Instead, they decode
# the 2-D graph first, then enumerate all possible 3-D forms and find the
# one with highest score.



msg = fn.copy_src(src="h", out="m")
in_dim=ATOM_FDIM

def collate(sample):
    graphs, labels = map(list,zip(*sample))
    batched_graph = dgl.batch(graphs)
    return batched_graph, torch.tensor(labels)

def reduce(nodes):
    # summazation by avarage is different part
    accum = torch.mean(nodes.mailbox['m'], 1)
    return {'h': accum}

class NodeApplyModule(nn.Module):
    def __init__(self, in_feats, out_feats, activation):
        super(NodeApplyModule, self).__init__()
        self.linear = nn.Linear(in_feats, out_feats)
        self.activation = activation
    
    def forward(self, node):
        h = self.linear(node.data['h'])
        h = self.activation(h)
        return {'h': h}


#fix bug
#https://discuss.dgl.ai/t/cant-run-gcn-chemo-py-an-example-on-github/1280/2
class GCN(nn.Module):
    def __init__(self, in_feats, out_feats, activation):
        super(GCN, self).__init__()
        self.linear = nn.Linear(in_feats, out_feats)
        self.activation = activation
    
    def forward(self, g, feature):
        g.ndata['h'] = feature
        g.update_all(msg, reduce)
        h = self.linear(g.ndata['h'])
        if self.activation is not None:
            h = self.activation(h)
        return h    
    
    
class Regressor(nn.Module):
    def __init__(self, in_dim, hidden_dim1,hidden_dim2, out_size,hidden_mode=False):
        super(Regressor, self).__init__()
        self.layers = nn.ModuleList([GCN(in_dim, hidden_dim1, F.relu),
                                    GCN(hidden_dim1, hidden_dim1, F.relu)])
        
        #two layer configuration
        self.regress1 = nn.Linear(hidden_dim1, hidden_dim2)
        self.regress2 = nn.Linear(hidden_dim2, out_size)
        
        self.hidden_mode=hidden_mode
    def forward(self, g):
        h = g.ndata['h']
        for conv in self.layers:
            h = conv(g, h)
        g.ndata['h'] = h
        hg = dgl.mean_nodes(g, 'h')
        hidden=self.regress1(hg)
        
        #output hidden layer
        if self.hidden_mode:
            return hidden
        
        return self.regress2(hidden)