import math
import matplotlib.pyplot as plt
def nHexagonFromCircle(n=1,R=1):
    m=6*n
    l1=[(0,R)]
    l2=[]
    l3=[(0,-R)]
    l4=[]
    p = (m//4 -1)if (m//4 == m/4) else m//4 
    for i in range(1,p+1):
        temp1=R*math.sin(2*math.pi*i/m)
        temp2=R*math.cos(2*math.pi*i/m)
        l1.append((temp1,temp2))
        l2.insert(0,(temp1,-temp2))
        l3.append((-temp1,-temp2))
        l4.insert(0,(-temp1,temp2))
    if (m//4 == m/4):
        l1.append((R,0))
        l3.append((-R,0))
    return l1+l2+l3+l4
def mesh_points(n=1,R=1):
    result=[]
    for i in range(1,n+1):
        result=result+nHexagonFromCircle(i,i*R/n)
    result.insert(0,(0,0))
    return result
def mesh_numRec(n):
    return 6*(n**2)
class cirle_mesh:
    def __init__(self,n=1,R=1):
        self.size=n
        self.radius=R
        self.points=mesh_points(n,R)
    def draw_mesh(self):
        x=[a[0] for a in self.points]
        y=[a[1] for a in self.points]
        for i in range(1,self.size+1): #ve cac da giac
            for j in range(1,i*6):
                id1=(i-1)*i*6//2+j
                plt.plot(x[id1:(id1+2)], y[id1:(id1+2)], 'b-')
            plt.plot([x[(i-1)*i*6//2+1], x[(i+1)*i*6//2]],[y[(i-1)*i*6//2+1], y[(i+1)*i*6//2]], 'b-')
        for i in range(2,self.size+1): # noi cac da giac
            id_i1=(i-1)*i*6//2+1
            id_i2=(i-2)*(i-1)*6//2+1
            for j in range(0,6):
                id_j1= j*i
                id_j2= id_j1-j
                for k in range(0,i):
                    id1= id_i1 + id_j1 + k
                    id2= id_i2 + id_j2 + k
                    if not((j == 5)and(k == (i-1))):
                        plt.plot([x[id1], x[id2]],[y[id1], y[id2]], 'b-')
                    else:
                        plt.plot([x[id1], x[(i-2)*(i-1)*6//2+1]],[y[id1], y[(i-2)*(i-1)*6//2+1]], 'b-')
                for k in range(1,i):
                    id1=(i-1)*i*6//2+1+j*i+k
                    id2=(i-2)*(i-1)*6//2+1+j*(i-1)+k-1
                    plt.plot([x[id1], x[id2]],[y[id1], y[id2]], 'b-')
        for i in range(1,7):
            plt.plot([x[0],x[i]],[y[0],y[i]], 'b-')
        plt.axis('equal')
        plt.show()
    def num_Rec(self):
        return mesh_numRec(self.size)
    def index_dinhke(self,index):
        if index > (self.count_iner_points()+self.size*6):
            print("self.index_dinhke(): index out of range")
            return []
        if index==0:
            return [t for t in range(1,7)]
        if index <= 6:
            l1=[0]
            if self.size == 1:
                l3=[]
            else:
                if index == 1:
                    l3=[7,8,18]
                else:
                    id1=6+(index-1)*2
                    l3=[id1,id1+1,id1+2]
            if index == 1:
                l2= [2,6]
            elif index == 6:
                l2= [1,5]
            else:
                l2= [index-1,index+1]
        else:
            i=1 
            while index > (i+1)*(i+2)*6//2:
                i=i+1
            id1= index- (i+1)*i*6//2 # i khong doi
            id2= id1-(id1//(i+1))*(i+1)
            if not(id2 == 1):
                id0=(i-1)*i*6//2
                if id2 !=0:
                    id3=id0 + (id1//(i+1))*i+id2-1
                    l1=[id3,id3+1]
                else:
                    id3=id0 + (id1//(i+1))*i+id2
                    if (id3+1) != (id3 + 6*i +1):
                        l1=[id3,id3+1]
                    else:
                        l1=[id0+1,id3]
                if i == (self.size-1):
                    l3=[]
                else:
                    id0=(i+1)*(i+2)*6//2
                    id3=id0 + (id1//(i+1))*(i+2)+id2
                    if id2 != 0:
                        l3=[id3, id3+1]
                    else:
                        l3=[id3-1,id3]
            else:
                if (id1//(i+1)) == 0:
                    id0=(i-1)*i*6//2
                    l1=[id0+1]
                    id3=id0+6*(2*i+1)
                    l3=[id3+1,id3+2,id3+6*(i+2)]
                else:
                    id0=(i-1)*i*6//2
                    l1=[id0+(id1//(i+1))*i+1]
                    id3= id0 + 6*(2*i+1) + (id1//(i+1))*(i+2) 
                    l3=[id3,id3+1,id3+2]
            if id1 == 1:
                l2= [index+1, (i+1)*(i+2)*6//2]
            elif id1 == 6*(i+1):
                l2= [index-6*(i+1)+1,index-1]
            else:
                l2= [index-1,index+1]
        return l1+l2+l3
    def cac_dinhke(self, index):
        return [self.points[i] for i in self.index_dinhke(index)]
    def all_triangle(self):
        l1=[(0,i,i+1) for i in range(1,6)]
        l1.append((0,1,6))
        if self.size == 1:
            return l1
        for i in range(2,self.size+1):
            id1=(i-1)*i*6//2
            for j in range(0,6):
                id2=i*j
                for k in range(1,i):
                    id3= id1 + id2 + k
                    l1.append(((i-2)*(i-1)*6//2+j*(i-1)+k,id3,id3+1))
                for k in range(2,i+1):
                    id3= (i-2)*(i-1)*6//2+j*(i-1)+k-1
                    if (id3+1) != (id1+1):
                        l1.append((id3,id3+1, id1+id2+k))
                    else:
                        l1.append(((i-2)*(i-1)*6//2+1,id1,id1+id2+k))
                if j != 5:
                    id3=id1+id2+i
                    l1.append(((i-2)*(i-1)*6//2+j*(i-1)+i,id3,id3+1))
                else:
                    l1.append(((i-2)*(i-1)*6//2+1,id1+1,id1+id2+i))
        return l1
    def count_iner_points(self):
        return 1+self.size*(self.size-1)*6//2

if __name__ == "__main__":
    a=cirle_mesh(10)
    #print(len(a.points)-6*a.size,a.count_iner_points())
    #index=2500
    #A=a.cac_dinhke(index)
    #print(a.num_Rec())
    #A=a.all_triangle()
    #print([a for a in A if not(a[0]<a[1]<a[2])])
    #plt.plot(a.points[index][0],a.points[index][1], 'ko')
    #x=[a[0] for a in A]
    #y=[a[1] for a in A]
    #plt.plot(x,y, 'ro')
    #plt.axis('equal')
    a.draw_mesh()
    #print(mesh_numRec(100))
    pass
