#object-oriented implementation

class node:
    def __init__(self, name, connections, branch_lengths):
        self.name = name
        self.connections = connections
        self.branch_lengths = branch_lengths
        self.num_con = len(self.connections)
        self.cbpairs = list(zip(connections,branch_lengths))
    def add_connection(self, connection, branch_length):      
        self.connections.append(connection)
        self.branch_lengths.append(branch_length)
        self.num_con = len(self.connections)
    def remove_connection(self,parent):
        self.cbpairs = [(x,y) for x,y in zip(self.connections,self.branch_lengths) if x!=parent]
        if len(self.connections) > 1:
            self.connections = list(list(zip(*self.cbpairs))[0])
            self.branch_lengths = list(list(zip(*self.cbpairs))[1])
        else:
            self.connections = []
            self.branch_length = []
    def set_parent(self,parent):
        self.parent = parent
    def expand(self):
        tmp1 = [a if isinstance(a,str) else 'intnode'+'%05.f'%a for a in self.connections]
        tmp2 = [str(a)+':'+'%g'%round(b,8) for a,b in zip(tmp1,self.branch_lengths)]
        return '('+','.join(tmp2)+')'
    def expand_rev(self):
        tmp1 = [a if isinstance(a,str) else 'intnode'+'%05.f'%a for a in self.connections]
        tmp2 = [str(a)+':'+'%g'%round(b,8) for a,b in zip(tmp1,self.branch_lengths)]
        return '('+','.join(tmp2[::-1])+')'

class newick:
    def __init__(self,text):
        self.text = text
        self.leaves = [x.split(':')[0].replace('(','') for x in text.split(',')]
        self.leaveswb = [x.split(')')[0].replace('(','') for x in text.split(',')]
        #remoe node labels
    def remove_node_labels(self):
        new = self.text
        c=0
        for i in self.text:
            if i == ')':
                c+=1
                new = ')'.join(new.split(')')[:c]) + ')' + ':' + ':'.join(')'.join(new.split(')')[c:]).split(':')[1:])     
        new = new[:-1]
        self.text = new

class tree:
    def __init__(self, nodes):
        self.nodes = nodes

    def add_node(self,node):
        self.nodes.append(node)

    def remove_node(self,name):
        self.nodes = [x for x in self.nodes if x.name!=name]

    def get_node(self, name):
        for i in self.nodes:
            if i.name == name:
                return i
        return None

    def remove_biconnection(self,name1,name2):
        self.get_node(name1).remove_connection(name2)
        self.get_node(name2).remove_connection(name1)
        
    def num_tips(self):
        c=0
        for i in self.nodes:
            if len(i.connections) == 1:
                c+=1
        return c

    def num_internal_nodes(self):
        return self.num_nodes() - self.num_tips()
        
    def num_nodes(self):
        return (len(self.nodes))
    
    def listnodes(self):
        return [x.name for x in self.nodes]
    def listtips(self):
        return [i.name for i in self.nodes if len(i.connections) == 1]
    
    def unresolved_nodes(self):
        for i in self.nodes:
            if len(i.connections) != 1 and len(i.connections) != 3: 
                print(i.name)
                print(i.connections)
            
    def root_at_tip(self,tip):
        assert len(self.get_node(tip).connections) == 1
        self.root = self.get_node(tip).connections[0]
    def root_at_node(self,nodename):
        self.root = nodename
        pass

    def export_nw(self,nt,parent):
        if nt == '':
            #initiate from root
            nt = self.get_node(self.root).expand()
            parent = self.root
        #check for internal nodes in nt
        for i in newick(nt).leaves:
            if 'intnode' in i:
                node_label = int(i.replace('intnode',''))
                self.remove_biconnection(parent,node_label)
                nt = nt.replace(i,self.get_node(node_label).expand())
                for i2 in self.get_node(node_label).connections: #remove all connections to and from expanded node
                    if not isinstance(i2,str):
                        self.remove_biconnection(node_label,i2)
                #print(nt,i,node_label,parent)
                #return self.export_nw(nt,node_label)
        for i in newick(nt).leaves:
            if 'intnode' in i:
                node_label = int(i.replace('intnode',''))
                #print(nt,i,node_label,parent)
                return self.export_nw(nt,node_label)
        return nt