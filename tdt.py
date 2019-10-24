#!/usr/bin/env python
import sys
import argparse
import logging
import time
import config
import netcon

# Global variables
TENSOR_NAMES = []
TENSOR_MATH_NAMES = []
BOND_NAMES = []
BOND_DIMS = []
VECTORS = []
FINAL_ORDER = None

class Tensor:
    def __init__(self, name = None, bonds = []):
        if name == None:
            self.name = []
        elif isinstance(name, list):
            self.name = name[:]
        else:
            self.name = [name]
        self.bonds = bonds[:]

    def __repr__(self):
        return "Tensor(" + str(self.name) + ", " + str(self.bonds) + ")"

    def __str__(self):
        return str(self.name) + ", " + str(self.bonds)


class Bond:
    def __init__(self, t0 = -1, t1 = -1):
        self.t0 = t0
        self.t1 = t1

    def __str__(self):
        return "({0}, {1})".format(self.t0, self.t1)

    def isFree(self):
        return (self.t0 < 0 or self.t1 < 0)

    def connect(self, tensor_index):
        assert self.isFree(), "edge already connected to two tensors"
        if self.t0 < 0:
            self.t0 = tensor_index
        else:
            assert not self.t0 == tensor_index, "edge connects to the same tensor"
            self.t1 = tensor_index

    def has(self, tensor_index):
        return (self.t0 == tensor_index or self.t1 == tensor_index)


class TensorNetwork:
    def __init__(self):
        self.tensors = []
        self.bonds = []
        self.total_memory = 0.0
        self.max_memory = 0.0
        self.cpu_cost = 0.0

    def __str__(self):
        s = ""
        for i, t in enumerate(self.tensors):
            s += "tensor {0} : {1}\n".format(i, t)
        for i, b in enumerate(self.bonds):
            s += "bond {0} : {1}, {2} {3}\n".format(i, BOND_NAMES[i], b, BOND_DIMS[i])
        s += "memory : {0}\n".format(self.total_memory)
        s += "cpu : {0}\n".format(self.cpu_cost)
        return s


    def clone(self):
        tn = TensorNetwork()
        tn.total_memory = self.total_memory
        tn.max_memory = self.max_memory
        tn.cpu_cost = self.cpu_cost
        tn.bonds = [Bond(b.t0, b.t1) for b in self.bonds]
        tn.tensors = [Tensor(t.name, t.bonds) for t in self.tensors]
        return tn


    def output_log(self, prefix=""):
        if not prefix == "": prefix += " "
        for i, t in enumerate(self.tensors):
            logging.info(prefix + "tensor{0} : {1} {2}".format(i, TENSOR_NAMES[i], t.bonds))
        for i,b in enumerate(self.bonds):
            logging.info(prefix + "bond{0} : {1} {2} {3}".format(i, BOND_NAMES[i], b, BOND_DIMS[i]))


    def add_tensor(self, t_name, b_names):
        t_index = len(self.tensors)
        b_indexs = []
        for b in b_names:
            if b not in BOND_NAMES:
                self.bonds.append(Bond())
                BOND_NAMES.append(b)
                BOND_DIMS.append(config.DEFAULT_BOND_DIM)

            i = BOND_NAMES.index(b)
            self.bonds[i].connect(t_index)
            b_indexs.append(i)

        TENSOR_NAMES.append(t_name)
        self.tensors.append(Tensor(t_index, b_indexs))


    def find_bonds(self, tensor_a, tensor_b):
        bonds_a = self.tensors[tensor_a].bonds
        bonds_b = self.tensors[tensor_b].bonds
        contract = [b for b in bonds_a if b in bonds_b]
        replaced_a = [b for b in bonds_a if b not in bonds_b]
        replaced_b = [b for b in bonds_b if b not in bonds_a]
        return contract, replaced_a, replaced_b


    def contract(self, t0, t1, bc, br0, br1):
        tn = self.clone()

        # create the contracted tensor
        t_new = tn.tensors[t0]
        ## change names of tensors using Reverse Polish Notation
        t_new.name = self.tensors[t0].name+self.tensors[t1].name+[-1]
        ## remove contracted bonds
        for b in bc: t_new.bonds.remove(b)
        ## add bonds from deleted tensor
        for b in br1: t_new.bonds.append(b)

        # clear the removed tensor
        tn.tensors[t1] = Tensor()

        # update bonds
        bonds = tn.bonds
        ## remove contracted bonds from the bond list
        for b in bc: bonds[b].t0 = bonds[b].t1 = -1
        ## change bond connections
        old_idx = t1
        new_idx = t0
        for b in br1:
            if bonds[b].t0 == old_idx: bonds[b].t0 = new_idx
            elif bonds[b].t1 == old_idx: bonds[b].t1 = new_idx

        return tn


    """Generate mathematical formula from Reverse Polish Notation"""
    stack = []
    for c in rpn:
        if c == -1:
            t1 = stack.pop()
            t0 = stack.pop()
            new_name = "(" + t0 + "*" + t1 + ")"
            stack.append(new_name)

        else:
            stack.append(TENSOR_MATH_NAMES[c])
    return stack[0]


#def get_connect_list_in(tn_orig, rpn, connect_list_in):
def get_connect_list_in(tn_orig, rpn):
    """Generate tensordot script from Reverse Polish Notation"""
    tn = tn_orig.clone()
    index = []
    name = []
    i = 1
    for c in rpn:
        #print(index)
        if c == -1:
            index1 = index.pop()
            index0 = index.pop()
            #print("index0", index0, "index1", index1)
            findices = sum([tn.tensors[index0].bonds,tn.tensors[index1].bonds], [])#flatten
            print([x for x in set(findices) if findices.count(x) > 1])
            for d in [x for x in set(findices) if findices.count(x) > 1]:
                for j in range(lenTN):
                    for k in range(len(connect_list_in[j])):
                        if connect_list_in[j][k] == BOND_NAMES[d]:
                            connect_list_in[j][k] = i
                i += 1
            #print(index)
            undep = [x for x in set(findices) if findices.count(x) == 1]
            tn.tensors[index0].bonds = undep
            tn.tensors[index1].bonds = []
            #t0 = tn.tensors[index0]
            #t1 = tn.tensors[index1]
            #bc, br0, br1 = tn.find_bonds(index0, index1)
            #axes0 = [t0.bonds.index(b) for b in bc]
            #axes1 = [t1.bonds.index(b) for b in bc]
            #tn = tn.contract(index0, index1, bc, br0, br1)
            #print(undep)
            #print(tn)
            #exit()
            index.append(index0)
        else:
            index.append(c)
    return
    #return connect_list_in


def read_file(infile, tn):
    """Read input file"""
    global FINAL_ORDER

    for line in infile:
        data = line.split()
        if data == []: continue

        command = data[0].lower()
        if command == "style":
            set_style(data[1].lower())
        elif command == "numpy":
            config.NUMPY = data[1]
        elif command == "indent":
            config.INDENT = " " * int(data[1])
        elif command == "default_dimension":
            # Should be set the top of input file.
            config.DEFAULT_BOND_DIM = int(data[1])
        elif command == "debug" or command == "verbose":
            config.LOGGING_LEVEL = logging.DEBUG
        elif command == "tensor":
            tn.add_tensor(data[1], data[2:])
        elif command == "bond":
            for b in data[1:-1]: set_bond_dim(b, int(data[-1]))
        elif command == "bond_dim":
            for b in data[2:]: set_bond_dim(b, int(data[1]))
        elif command == "order":
            FINAL_ORDER = data[1:]
        elif command == "vector":
            VECTORS.append((data[1], data[2]))
    infile.close()


def check_bond_order(tn):
    return FINAL_ORDER == None or \
        frozenset(FINAL_ORDER) == frozenset(BOND_NAMES[i] for i, b in enumerate(tn.bonds) if b.isFree())


def check_vector():
    for v in VECTORS:
        if v[1] not in BOND_NAMES: return False
    return True


def parse_args():
    parser = argparse.ArgumentParser(description = "Code generator for tensor contraction")
    parser.add_argument('-s', metavar = 'style', dest = 'style',
                        type = str, default = None,
                        choices = ['numpy', 'mptensor', 'julia'],
                        help = 'set output style ("numpy" or "mptensor" or "julia")')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose',
                        help = 'verbose mode')
    parser.add_argument('-o', metavar = 'outfile', dest = 'outfile',
                        type = argparse.FileType('w'), default = sys.stdout,
                        help = 'write the result to outfile')
    parser.add_argument('infile',
                        type = argparse.FileType('r'),
                        help = 'tensor-network definition file')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    tn = TensorNetwork()

    # Read input file
    read_file(args.infile, tn)

    # Overwrite by command-line option
    set_style(args.style)
    if args.verbose:
        config.LOGGING_LEVEL = logging.DEBUG

    assert len(tn.tensors) > 0, "No tensor."
    assert len(tn.bonds) > 0, "No bond."
    assert check_bond_order(tn), "Final bond order is invalid."
    assert check_vector(), "Vectors will be put on non-existent bond."
    logging.basicConfig(format = "%(levelname)s:%(message)s", level = config.LOGGING_LEVEL)
    tn.output_log("input")
    rpn, cpu = netcon.NetconOptimizer(tn.tensors, BOND_DIMS).optimize()
    
    #set tensor_list
    tensor_list = "Any["
    lenTN = len(TENSOR_NAMES)
    for i in range(lenTN):
        tensor_list += TENSOR_NAMES[i]
        if i != lenTN - 1: 
            tensor_list += ", "
    tensor_list += "]"

    connect_list_in = [[BOND_NAMES[i] for i in tn.tensors[j].bonds] for j in range(len(tn.tensors))]
    
    #set FINAL_ORDER if it is None
    if FINAL_ORDER == None:
        fcli = sum(connect_list_in, []) #list flatten
        FINAL_ORDER = [x for x in set(fcli) if fcli.count(x) == 1]

    #set output minus index
    for i in range(len(FINAL_ORDER)):    
        for j in range(lenTN):
            for k in range(len(connect_list_in[j])):
                if connect_list_in[j][k] == FINAL_ORDER[i]:
                    connect_list_in[j][k] = -i-1
            else:
                continue
            break
        else:
            continue
        break       

    get_connect_list_in(tn, rpn)
    print(connect_list_in)
