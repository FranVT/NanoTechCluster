"""
    Script to create all pairwise tables interactions
"""

include("functions.jl")

## 3 body interaction

# Parameters
N = 2^5;
M = 2*N*N*N;

eps_ij = 10.0;
eps_ik = 10.0;
eps_jk = 10.0;
sig = 0.4;
rmin = sig-sig/10;
rmax = 1.499*sig;
thi = 180/(4*N)
thf = 180 - thi;

# Create the domains of evaluation according filename nessetities
th_dom = range(thi,thf,2*N);
r_dom = range(rmin,rmax,N);

doms = reduce(vcat,Iterators.product(r_dom,r_dom,th_dom));

# Create tuples with the information
docs =  map(eachindex(doms)) do s
            (
                 s,
                 doms[s]...,
                 -Forceij(eps_ij,eps_ik,eps_jk,sig,doms[s][1],doms[s][2]),
                 -Forceik(eps_ij,eps_ik,eps_jk,sig,doms[s][1],doms[s][2]),
                 Forceij(eps_ij,eps_ik,eps_jk,sig,doms[s][1],doms[s][2]),
                 0.0,
                 Forceik(eps_ij,eps_ik,eps_jk,sig,doms[s][1],doms[s][2]),
                 0.0,
                 SwapU(1.0,eps_ij,eps_ik,eps_jk,sig,doms[s][1],doms[s][2])
            )
        end

# Create the 3 body table file
createTable(N,rmin,rmax,docs,"swapMechTab.table");

# Create the swap.3b file
aux = vec(Iterators.product(Iterators.repeated(("PA ","PB "),3)...)|>collect);
parameters = string("0.6 swapMechTab.table SEC1 linear ",N,"\n");
info = join.(aux).*parameters;

# Start to write the data file
touch("swapMech.3b"); # Create the file

# Edit the file
open("swapMech.3b","w") do f
    map(s->write(f,s),info)
end

## Pairwise patch-patch interaction

# Parameters
N = 1000000; #2^20;
sig = 0.4;
eps = 1.0;
rmin = 0.000001;
rmax = 5*sig;
r_dom = range(rmin,rmax,length=N);

# Create the table
info = map(s->(s,r_dom[s],Upatch(eps,sig,r_dom[s]),Fpatch(eps,sig,r_dom[s])),eachindex(r_dom));

# Start to write the data file
touch("patchTab.table"); # Create the file

# Edit the file
open("patchTab.table","w") do f
    write(f,"DATE: 2024-07-08 UNITS: lj CONTRIBUTOR: Fco.\n\n\n")
    write(f,"POT\n")
    write(f,string("N ",N,"\n\n"))
    map(t->write(f,rstrip(join(map(s->s*" ",string.(info[t]))))*"\n" ),eachindex(info))
end





