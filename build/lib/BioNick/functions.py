#import random
import pandas as pd
import matplotlib.pyplot as plt

#test data
#bt  = "(t6:0.3806255851,(t7:0.5440872659,(t5:0.6203179485,t1:0.7089423968):0.3421749519):0.2618428892,((t4:0.5886056426,t3:0.8614832051):0.09094934678,(t8:0.6766210056,t2:0.1412485428):0.4451906278):0.2356684168)"
#ot = '(((((((((A_O.sativa:0.1,A_O.glaberrima:0.1):0.1,(A_O.barthii:0.1,A_O.glumipatula:0.1):0.1):0.1,(A_O.meridionalis:0.1,A_O.nivara:0.1,A_O.rufipogon:0.1):0.1):0.1,B_O.punctata:0.1):0.1,((C_O.officinalis:0.1,C_O.alta:0.1):0.1,D_O.alta:0.1):0.1):0.1,E_O.australiensis:0.1):0.1,F_O.brachyantha:0.1):0.1,(K_O.coarctata:0.1,L_O.coarctata:0.1):0.1):0.1,OG_L.perrieri:0.1)'
#ot2 = '((((((((A_O.sativa:0.1,A_O.glaberrima:0.1):0.1,(A_O.barthii:0.1,A_O.glumipatula:0.1):0.1):0.1,(A_O.meridionalis:0.1,A_O.nivara:0.1,A_O.rufipogon:0.1):0.1):0.1,B_O.punctata:0.1):0.1,((C_O.officinalis:0.1,C_O.alta:0.1):0.1,D_O.alta:0.1):0.1):0.1,E_O.australiensis:0.1):0.1,F_O.brachyantha:0.1):0.1,(K_O.coarctata:0.1,L_O.coarctata:0.1):0.1,OG_L.perrieri:0.1)'


##########################################################################


#extract leaves
def leaves(nw):
    return [x.split(':')[0].replace('(','') for x in nw.split(',')]

#extract leaves with branch length
def leaves_wb(nw):
    return [x.split(')')[0].replace('(','') for x in nw.split(',')]

#make sure labels are unique. Append numbers if they are not
#assert len(leaves(bt)) == len(set(leaves(bt)))

#label leaves
def lale(nw):
    t = {}
    for i in leaves(nw):
        t[i] = leaves(nw).index(i)
    return t

#convert nw to tsv
def recur_nw_pd(bt,n,t): #,n = len(leaves(bt)),t = []): #default arguments cause trouble
    c = 0
    for i in bt:
        if i == '(':
            c+=1
        if i == ')':
            nt = bt.replace('('+bt.split('(')[c].split(')')[0]+')', '__'+str(n))
            for i2 in bt.split('(')[c].split(')')[0].split(','):
                t.append([n,i2.split(':')[0], float(i2.split(':')[1])])
            n+=1
            c-=1
            return recur_nw_pd(nt,n,t)
    return bt,t


#convert nw to tsv
def nw_pd(tree):
    a,b = recur_nw_pd(tree,len(leaves(tree)),[])
    return b

#replace leaves with labels
def encode_leaves(bt,ab):
    #replace leaves with labels
    tb = ab
    ll = lale(bt)
    n = 0
    for a,b,c in tb:
        if b in ll.keys():
            tb[n][1] = ll[b]
        else:
            tb[n][1] = int(b.replace('__',''))
        n+=1
    return tb


#swap root node
def swap_root(tb,x):
    df = pd.DataFrame(tb)
    df2 = df.copy()
    while (x in df[1].unique()):
        #swap
        tmp = int(df.loc[df[1] == x, 0].iloc[0])
        df2.loc[df[1] == x, 0] = int(df2.loc[df[1] == x, 1].iloc[0])
        df2.loc[df[1] == x, 1] = tmp

        #print(tmp,x)
        x = tmp
    return df2


##convert tsv to newick.

#add 10 trailing decimals to the numbers.
def trail(df):
    t = []
    for a,b,c in df.values:
        if isinstance(a,int):
            na = '%016.10f' % a
        else:
            na = a
        if isinstance(b,int):
            nb = '%016.10f' % b
        else:
            nb = b
        t.append((na,nb,c))   
    return pd.DataFrame(t)

#reassign leaf names
def reasign(df2,bt):
    ll = lale(bt)
    inv = {v: k for k, v in ll.items()} #invert keys
    t = []
    for i in df2[1]:
        if i in inv.keys():
            t.append(inv[i])
        else:
            t.append(i)
    df2[1] = t

    df3 = df2.copy()
    
    return df3

def expand_node(df,node):
    tmp = df[df[0] == '%016.10f' % float(node)] #dependency on pandas.
    
    expanded_node = ''
    for a,b,c in tmp.values: #tmp.values when depends on pandas
        expanded_node = expanded_node+','+str(b)+':'+str(c)
    expanded_node = expanded_node[1:]
    return '(' + expanded_node + ')'

def recur_pd_nw(nt,df):
    #lvs = set(list(zip(*tb))[0]) #internal nodes
    int_nodes = df[0].astype(str).unique()
    for i in leaves(nt):
        if i in int_nodes:
            return recur_pd_nw(nt.replace(i,expand_node(df,i)),df)
    return nt    

#root at taxon
def root_at(tree,taxon):
    a,b = recur_nw_pd(tree,len(leaves(tree)),[])
    d = pd.DataFrame(b)
    i = d.loc[d[1] == taxon,0].iloc[0] #new root    
    tb = encode_leaves(tree,b)
    return recur_pd_nw('%016.10f' % i, trail(reasign(swap_root(tb,i),tree)))

#root at node
def root_at_node(tree,i):
    a,b = recur_nw_pd(tree,len(leaves(tree)),[])    
    tb = encode_leaves(tree,b)
    return recur_pd_nw('%016.10f' % i, trail(reasign(swap_root(tb,i),tree)))

#flip label order
def expand_node_flip(df,node):
    tmp = df[df[0] == '%016.10f'% float(node)] #dependency on pandas.
    expanded_node = ''
    for a,b,c in tmp.values: #tmp.values when depends on pandas
        expanded_node = str(b)+':'+str(c)+','+expanded_node
    expanded_node = expanded_node[:-1]
    return '(' + expanded_node + ')'
def recur_pd_nw_flip(nt,df):
    int_nodes = df[0].astype(str).unique()
    for i in leaves(nt):
        if i in int_nodes:
            return recur_pd_nw_flip(nt.replace(i,expand_node_flip(df,i)),df)
    return nt    

def flip_all_edges(tree):
    root_node,b = (recur_nw_pd(tree,len(leaves(tree)),[]))
    root_node = int(root_node.replace('__',''))
    tb = encode_leaves(tree,b)
    new = recur_pd_nw_flip('%016.10f' % root_node, trail(reasign(pd.DataFrame(tb),tree))) 
    return new

#flip specified node edges
def recur_pd_nw_flip_at_node(nt,df,node):
    int_nodes = df[0].astype(str).unique()
    for i in leaves(nt):
        if i in int_nodes:
            if i == str('%016.10f' % node):
                return recur_pd_nw_flip_at_node(nt.replace(i,expand_node_flip(df,i)),df,node)
            else:
                return recur_pd_nw_flip_at_node(nt.replace(i,expand_node(df,i)),df,node)  
    return nt   
def flip_leaves_at_node(tree,node):
    root_node,b = (recur_nw_pd(tree,len(leaves(tree)),[]))
    root_node = int(root_node.replace('__',''))
    tb = encode_leaves(tree,b)
    new = recur_pd_nw_flip_at_node('%016.10f' % root_node, trail(reasign(pd.DataFrame(tb),tree)),node)
    return new


#travel through pandas 
## A travelling tree string will still not be identical. The float conversion truncates trailing decimal zeros.
def travel(tree):
    root_node,b = (recur_nw_pd(tree,len(leaves(tree)),[]))
    root_node = int(root_node.replace('__',''))
    tb = encode_leaves(tree,b)
    new = recur_pd_nw('%016.10f' % root_node, trail(reasign(pd.DataFrame(tb),tree)))
    return new


#export trees rooted at every internal node
def all_trees(bt):
    a,b = recur_nw_pd(bt,len(leaves(bt)),[])
    x = set(list(zip(*b))[0])
    
    tb = encode_leaves(bt,b)
    
    t = []
    for i in x:
        t.append(recur_pd_nw('%016.10f' % i, trail(reasign(swap_root(tb,i),bt))))
        print(recur_pd_nw('%016.10f' % i, trail(reasign(swap_root(tb,i),bt))))
    return t

#export all node descendents to a dictionary 
def nodes_w_all_descendants(tree):
    t2 = {}
    n = 0
    for i in leaves(tree):
        t2[i] = [i]
        n+=1
    a,b,c = recur_nw_pd_an(tree,len(leaves(tree)),[],t2)
    return a,b,c
    
def recur_nw_pd_an(bt,n,t,t2):
    c = 0
    t2['__'+str(n)] = []
    for i in bt:
        if i == '(':
            c+=1
        if i == ')':
            nt = bt.replace('('+bt.split('(')[c].split(')')[0]+')', '__'+str(n))
            for i2 in bt.split('(')[c].split(')')[0].split(','):
                t.append([n,i2.split(':')[0], float(i2.split(':')[1])])
                t2['__'+str(n)] = t2['__'+str(n)] + t2[i2.split(':')[0]]
            n+=1
            c-=1
            return recur_nw_pd_an(nt,n,t,t2)
    return bt,t,t2

def extract_subtree(tree, leaves_to_keep):
    t = tree
    for leaf in leaves(t):
        if leaf not in leaves_to_keep:
            t = remove_leaf(t,leaf)
    a,b = recur_nw_pd(t,len(leaves(t)),[])
    #prune root
    bd = pd.DataFrame(b)
    mi = bd[bd[0].duplicated(keep = False)][0].max()
    #mi = bd[bd[0] == mn].iloc[-1,0]
    b2 = []
    for l,r,bl in b:
        if l > mi:
            break
        b2.append([l,r,bl])
    bd = pd.DataFrame(b2)
    
    #if nodes don't need joining
    if bd[~bd[0].duplicated(keep = False)].shape[0] == 0:
        c = pd.DataFrame(encode_leaves(t,b2))
        c = trail(reasign(c,t))
        return recur_pd_nw('%016.10f' % mi,c)
    
    #join_singular_nodes
    carry_over = 0
    pl = bd[~bd[0].duplicated(keep = False)].values[0][0]
    pr = bd[~bd[0].duplicated(keep = False)].values[0][1]
    for l,r,bl in bd[~bd[0].duplicated(keep = False)].values[1:]:
        if pl == l-1:
            #print(l,r,bl)
            carry_over = carry_over + bl
        else:
            bd.loc[bd[1] == '__'+str(pl),2] = bd.loc[bd[1] == '__'+str(pl),2] + carry_over
            bd.loc[bd[1] == '__'+str(pl),1] = pr
            pr = r
        pl = l
        #pr = r
    bd.loc[bd[1] == '__'+str(pl),2] = bd.loc[bd[1] == '__'+str(pl),2] + carry_over
    bd.loc[bd[1] == '__'+str(pl),1] = pr
    
    b2 = bd[bd[0].duplicated(keep = False)].round(decimals=8).values
    
    c = pd.DataFrame(encode_leaves(t,b2))
    c = trail(reasign(c,t))
    
    return recur_pd_nw('%016.10f' % mi,c)

#singularize nodes in a tree. i.e., remove extra nodes.
def recur_sin(new):
    for leaf in leaves_wb(new):
        enclosed_leaf = '('+leaf+')'
        if enclosed_leaf in new:
            leaf_without_branch = leaf.split(':')[0]
            inner_branch = float(leaf.split(':')[1])
            outer_branch = float(new.split(enclosed_leaf)[1].split(',')[0].split(')')[0][1:])
            new_branch = round(inner_branch+outer_branch,9)
            enclosed_leaf_with_outer_branch = enclosed_leaf+':'+f"{outer_branch:.10f}"
            new = new.replace(enclosed_leaf_with_outer_branch, leaf_without_branch+':'+f"{new_branch:.10f}")
            return recur_sin(new)
    return new

#remove a leaf from a tree
def remove_leaf(tree,name):
    name = name+':'
    leading = tree.split(name)[0]
    if ',' in tree.split(name)[1]:
        comma = tree.split(name)[1].index(',')
    else:
        comma = 0
    close = tree.split(name)[1].index(')')
    if close < comma:
        trailing = ')'.join(tree.split(name)[1].split(')')[1:])
        new = leading[:-1] + ')' + trailing #:-1 removes last comma
    else:
        trailing = ','.join(tree.split(name)[1].split(',')[1:])
        new = leading + trailing
    
    new = recur_sin(new)
    
    if new[-1] == ',': #if terminal leaf removed
        new = new[:-1]+')'
        
    return new


#draw a cladogram
def draw_clad(tree):
    root_node,b = recur_nw_pd(tree,len(leaves(tree)),[])    
    m = pd.DataFrame(encode_leaves(tree,b))
    
    #vertical node locations
    t1 = dict([(i,i) for i in range(len(leaves(tree)))])
    for a,b,c in m.values:
        if a not in t1.keys():
            t1[a] = mean_d(m,a,t1)
    
    #horizontal node locations. end of edge. point at the right.      
    t2 = {int(root_node.replace('__','')):0}
    for a,b,c in m.values[::-1]:
        if b not in t2.keys():
            t2[b] = t2[a]+c
            
    #horizontal node locations. start of edge. point at the left.  
    t3 = {int(root_node.replace('__','')):0}
    for a,b,c in m.values[::-1]:
        t3[b] = t2[a]
    
    #vertical start and end.
    t4 = {}
    for a,b,c in m.values:
        if a not in t4.keys():
            t4[a] = ab(m,a,t1)
            
    n = pd.concat([pd.DataFrame(t1.values(),index = t1.keys()), 
                   pd.DataFrame(t2.values(),index = t2.keys()), 
                   pd.DataFrame(t3.values(),index = t3.keys())], axis = 1)    
    
    #nodes
    #plt.scatter(n.iloc[:,1],n.iloc[:,0])
    
    #horizontal lines
    for a,b,c in n.values:
        plt.plot([b,c],[a,a], color = 'darkblue')
    
    #vertical lines
    for k,v in t4.items():
        plt.plot([t2[k],t2[k]],v, color = 'darkblue')

    #dash to leaves
#     for a,b,c in n.values[:len(leaves(tree))]:
#         plt.plot([c,max(t2.values())],[a,a], linestyle = '--', color = 'darkblue')
            
#for vertical distancing of nodes. mean of first and last edge that go right from a node.
def mean_d(m,node,t1):
    sw = False
    for a,b,c in m.values:
        if (sw == False) and (a == node):
            sw = True
            s = t1[b]
        if sw and a != node:
            return (s+pb)/2
        pb = t1[b]
    return (s+pb)/2

#for vertical lines. first and last edge that go right from a node.
def ab(m,node,t1):
    sw = False
    for a,b,c in m.values:
        if (sw == False) and (a == node):
            sw = True
            s = t1[b]
        if sw and a != node:
            return [s,pb]
        pb = t1[b]
    return [s,pb]