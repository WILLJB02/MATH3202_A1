#Given Data Stubs
IDtoLVC = [
[14.4,55.2,62.1,57.0,42.4,80.7,63.7,8.8],
[46.6,5.9,58.1,78.2,77.7,55.3,14.4,63.6],
[66.7,68.0,8.9,24.0,39.6,29.2,55.9,61.5]
]

CCDPop = [4030,3355,3920,4271,4017,2382,3935,2316,4322,2151,3774,3699,4533,3734,2992,3624,4755,4045,4035,2238,4661,3026,3570,4090,2548]

CCDtoLVC = [
[19.9,0,0,0,0,0,0,18.0],
[7.7,0,0,0,0,0,0,16.9],
[18.4,23.7,0,0,0,0,0,0],
[0,14.8,0,0,0,0,33.5,0],
[0,14.7,0,0,0,0,20.6,0],
[0,0,0,0,28.7,0,0,15.4],
[8.8,0,0,0,0,0,0,15.1],
[27.7,23.5,0,0,0,0,25.9,0],
[0,13.2,0,0,0,0,18.5,0],
[0,25.9,0,0,0,0,14.4,0],
[0,0,0,30.7,14.2,0,0,19.7],
[30.8,0,0,0,20.4,0,0,19.3],
[0,0,21.2,0,0,0,35.2,0],
[0,0,29.2,0,0,27.2,17.9,0],
[0,0,0,0,0,27.7,23.5,0],
[0,0,0,11.4,9.5,0,0,0],
[0,0,24.5,5.9,17.7,0,0,0],
[0,0,10.8,18.8,29.7,0,0,0],
[0,0,20.9,0,0,12.5,31.2,0],
[0,0,0,0,0,14.0,35.2,0],
[0,0,0,18.8,27.5,0,0,0],
[0,0,0,9.4,25.2,0,0,0],
[0,0,16.8,28.5,0,30.0,0,0],
[0,0,24.0,0,0,21.9,0,0],
[0,0,36.6,0,0,14.9,0,0]
]

from gurobipy import *

#Creating Model
m = Model("VaccineDistribution")

#Sets
ID = ["ID-A","ID-B","ID-C"]
LVC = [("LVC" + (str(i))) for i in range(8)]
CCD = [str(i) for i in range(25)]
D = {(ID[IDindex], LVC[LVCindex]) : DistL 
     for IDindex, lenghtList in enumerate(IDtoLVC) 
     for LVCindex, DistL in enumerate(lenghtList)}
R = {(CCD[CCDindex], LVC[LVCindex]) : DriveL
     for CCDindex, lengthList in enumerate(CCDtoLVC)
     for LVCindex, DriveL in enumerate(lengthList)
     if DriveL != 0}           
T = [i for i in range(6)]

#Data
ImportC = {"ID-A": 177, "ID-B": 142, "ID-C": 155}
P = {CCD : CCDPop[index] for index, CCD in enumerate(CCD)}
DelC = 0.2
"DelL_d = D[d] for d element of D"
"o_d = d[0] for d element of D"
"v_d = d[1] for d element of D"
"DriveL = R[r] for r element of R"
DriveC = 1
"t_r = r[0] for r element of R"
"f_r = r[1] for r element of R"
Cap = 37000
TMax = 16000
WMax = 2200
DelayC = 10
VaxDiff = 0.1

#Variables
X = {(d,t): m.addVar() for d in D for t in T}
Y = {i: m.addVar() for i in ID}       
Z = {(r,t): m.addVar() for r in R for t in T}   
S = {(l,t): m.addVar() for l in LVC for t in T}
UNVAC = {(c,t): m.addVar() for c in CCD for t in T}

#Objective
m.setObjective(quicksum(ImportC[i]*Y[i] for i in ID) + 
               quicksum(DelC * D[d] * X[d,t] for d in D for t in T) + 
               quicksum(DriveC * R[r] * Z[r,t] for r in R for t in T) +
               quicksum(DelayC * UNVAC[c,t] for c in CCD for t in T), GRB.MINIMIZE)

#Constraints
Import_Cap =[]
for i in ID:
    c = m.addConstr(Y[i] <= Cap)
    Import_Cap.append(c)
    m.addConstr(quicksum(X[d,t] for d in D for t in T if d[0] == i) <= Y[i])
    
for l in LVC:
    for t in T:
        if t > 0:
            m.addConstr(S[l,t] == S[l,t-1] + 
                        quicksum(X[d,t] for d in D if d[1] == l) - 
                        quicksum(Z[r,t] for r in R if r[1] == l))
        else:
             m.addConstr(S[l,t] ==  
                        quicksum(X[d,t] for d in D if d[1] == l) - 
                        quicksum(Z[r,t] for r in R if r[1] == l))

Weekly_Vax_Cap =[]
for l in LVC:
        m.addConstr(quicksum(Z[r,t] for r in R for t in T if r[1] == l) <= TMax)
        for t in T:
            c = m.addConstr(quicksum(Z[r,t] for r in R if r[1] == l) <= WMax)
            Weekly_Vax_Cap.append(c)

for c in CCD:
        m.addConstr(quicksum(Z[r,t] for r in R for t in T if r[0] == c) == P[c])
        for t in T:
            if t > 0:
                m.addConstr(UNVAC[c,t] == UNVAC[c,t-1] - quicksum(Z[r,t] for r in R if r[0] == c))
            else:
                m.addConstr(UNVAC[c,t] == P[c] - quicksum(Z[r,t] for r in R if r[0] == c))

for c1 in CCD:
    for c2 in CCD:
        for t in T:
            m.addConstr(-VaxDiff <= UNVAC[c1,t] / P[c1] - UNVAC[c2,t] / P[c2])
            m.addConstr(VaxDiff >= UNVAC[c1,t] / P[c1] - UNVAC[c2,t] / P[c2])


#Optimizing
m.optimize()


#Showing Results
print("\n\n----------------------------------------------------")
print("Optimized Distribution Cost: " + str(m.objVal))
print("----------------------------------------------------")
print("Number of Vaccines imported to IDs:")
for i in ID:
    print(str(i) + ": " + str((Y[i].x)))
print("----------------------------------------------------")   
print("Vaccine Distribution from ID to LVC:")
print('{0:<20} {1:<10} {2}'.format("(ID, LVC)", "TotaL", "Weekly"))
for d in D:
    weekly_doses = [round((X[d,t].x),1) for t in T]
    if sum(weekly_doses) != 0:
        print('{0:<20} {1:<10} {2}'.format(str(d), 
                                           str(round(sum(weekly_doses),1)), 
                                           str(weekly_doses)))
print("----------------------------------------------------")   
print("Vaccine Adminstration:")
print('{0:<20} {1:<10} {2}'.format("LVC", "TotaL", "Weekly"))
for l in LVC:
    weekly_administer = []
    for t in T:
        vaccines = 0
        for r in R:
            if r[1] == l:
                vaccines = vaccines + Z[r,t].x
        weekly_administer.append(round(vaccines,1))
    if sum(weekly_administer) != 0:
        print('{0:<20} {1:<10} {2}'.format(str(l), 
                                           str(round(sum(weekly_administer),1)), 
                                           str(weekly_administer)))

print("----------------------------------------------------")   
print("Vaccine Take Up from CCD to LVC:")
print('{0:<20} {1:<10} {2}'.format("(ID, LVC)", "TotaL", "Weekly"))
for r in R:
    weekly_takeup = [round((Z[r,t].x),1) for t in T]
    if sum(weekly_takeup) != 0:
        print('{0:<20} {1:<10} {2}'.format(str(r), 
                                           str(round(sum(weekly_takeup),1)), 
                                           str(weekly_takeup)))
        
print("----------------------------------------------------")  
print("Percentage of Unvaccinated People:")
for c in CCD:
        print(str(c) + ": " + str([round(UNVAC[c,t].x/P[c]*100, 1) for t in T]))
print("----------------------------------------------------")     
print("Improving Weekly Vaccine Administration Cap:")
print('{0:<10} {1:<10} {2:<10} {3}'.format("LVC", "Decrease", "Lower", "Upper"))
Dual = []
lower_bound = []
upper_bound = []
index = 0
for c in Weekly_Vax_Cap:
    if index % 6 == 0: 
        dual = 0
        low = []
        u = []
    dual = dual + c.Pi
    low.append(c.SARHSLow)
    u.append(c.SARHSUp)  
    if index % 6 == 5:
        print('{0:<10} {1:<10} {2:<10} {3}'.format(LVC[int(index/6)], 
                                                   round(dual), 
                                                   round(max(low)), 
                                                   round(min(u))))
        Dual.append(dual)
        lower_bound.append(max(low))
        upper_bound.append(min(u)) 
    index = index + 1 
print('{0:<10} {1:<10} {2:<10} {3}'.format("Total", round(sum(Dual)), 
                                           str(round(max(lower_bound))), 
                                           str(round(min(upper_bound)))))
print("----------------------------------------------------")     
print("Improving Importation Cap:")   
Dual = []
lower_bound = []
upper_bound = []
for i, c in enumerate(Import_Cap):
    if i != 0:
        print('{0:<10} {1:<10} {2:<10} {3}'.format(ID[i], 
                                                   round(c.Pi), 
                                                   round(c.SARHSLow), 
                                                   round(c.SARHSUp)))
    else:
        print('{0:<10} {1:<10} {2:<10} {3}'.format(ID[i], 
                                                   round(c.Pi), 
                                                   round(c.SARHSLow), 
                                                   c.SARHSUp))
    Dual.append(c.Pi)
    lower_bound.append(c.SARHSLow)
    upper_bound.append(c.SARHSUp) 
print('{0:<10} {1:<10} {2:<10} {3}'.format("Total", 
                                           round(sum(Dual)), 
                                           str(round(max(lower_bound))), 
                                           str(round(min(upper_bound)))))




