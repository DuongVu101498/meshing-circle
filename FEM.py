import numpy as np
import math
import pylab as pl
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
from matplotlib import cm
import cmesh
import time
from scipy.integrate import dblquad
GAUSS =((1/3, 1/3, 9/40), ((6-math.sqrt(15))/21, (6-math.sqrt(15))/21, (155-math.sqrt(15))/1200),
((9+2*math.sqrt(15))/21, (6-math.sqrt(15))/21, (155-math.sqrt(15))/1200),
((6-math.sqrt(15))/21, (9+2*math.sqrt(15))/21, (155-math.sqrt(15))/1200),
((6+math.sqrt(15))/21, (6+math.sqrt(15))/21, (155+math.sqrt(15))/1200),
((6+math.sqrt(15))/21,(9-2*math.sqrt(15))/21,(155+math.sqrt(15))/1200),
((9-2*math.sqrt(15))/21, (6+math.sqrt(15))/21,(155+math.sqrt(15))/1200))

f1= lambda x,y: 8*(x+y)  
u1= lambda x,y:(x+y)*(1-x**2-y**2)
f2= lambda x,y: 4
u2= lambda x,y: 1-x**2-y**2
f3= lambda x,y: -0.1*(-22*x**3-6*x*y**2-14*y**2-2*x**2-10*x-24*y-30)
u3= lambda x,y: 0.1*(x**3+y**2+2*x+3*y+8)*(1-x**2-y**2)
f4= lambda x,y: -20*x*y*(y**2-x**2)*20
u4= lambda x,y: x*y*(x**2-y**2)*(1-x**2-y**2)*20
f5= lambda x,y: (1-x**2-y**2)**-1.5+(1-x**2-y**2)**-0.5
u5= lambda x,y: (1-x**2-y**2)**0.5 if 1-x**2-y**2>0 else 0
f6= lambda x,y: -0.1*(16*math.pi*math.cos(4*math.pi*(x**2+y**2-1))-64*math.pi**2*(x**2+y**2)*math.sin(4*math.pi*(x**2+y**2-1)))
u6= lambda x,y: 0.1*math.sin(4*math.pi*(x**2+y**2-1))
f7= lambda x,y: -4*(x**2+y**2-1)*math.exp(1-x**2-y**2)/(math.exp(1)-1)
u7= lambda x,y: (math.exp(1-x**2-y**2)-1)/(math.exp(1)-1)

_resolution=16
_f=f2
_u=u2
_toa_do = (0,0)
_Radius=1
_draw_nghiem= 1
_draw_u=0
class virtual_symmetric_matrix:
    def __init__(self, size):
        self.shape= (size, size)
        self.data= [[[i,0]] for i in range(0,size)]
    def set(self, i, j, value):
        if j<i or i> self.shape[0] or j>self.shape[0]:
            print("set method Error: i must <= j , or i, j out of range")
        else:
            if j not in [a[0] for a in self.data[i]]:
                self.data[i].append([j,value])
                self.data[j].append([i,value])
            else:
                for a in  self.data[i]:
                    if a[0] == j:
                        a[1]=value
                        break
                for a in  self.data[j]:
                    if a[0] == i:
                        a[1]=value
                        break
    def get(self, i, j):
        if i <= j:
            for a in  self.data[i]:
                    if a[0] == j:
                        return a[1]
            return 0
        else:
            return self.get(j,i)
def CGM(vs_matrix,b):
    start_time_CGM = time.time()
    x0=[0 for i in range(0,vs_matrix.shape[0])]
    r0=list(b)
    p0=list(b)
    k=0
    while k<= vs_matrix.shape[0]:
        tmp1=0 # (r_k)T.r_k
        tmp2=0 #  (p_k)T.A.p_k
        for i in range(0,vs_matrix.shape[0]):
            tmp1=tmp1+ r0[i]**2
        for i in range(0,vs_matrix.shape[0]):
            for t in vs_matrix.data[i]:
                if t[0] == i:
                    tmp2=tmp2+t[1]*p0[i]*p0[t[0]]
                elif t[0] > i:
                    tmp2=tmp2+2*t[1]*p0[i]*p0[t[0]]
        alpha = tmp1/tmp2 # alpha
        x1= [x0[i]+ alpha*p0[i] for i in range(0,vs_matrix.shape[0])] # x_k+1
        r1=[0 for i in range(0,vs_matrix.shape[0])]                   # r_k+1
        for i in range(0,vs_matrix.shape[0]):
            tmp3=0
            for t in vs_matrix.data[i]:
                tmp3=tmp3+t[1]*p0[t[0]]
            r1[i]=r0[i]-alpha*tmp3
        ### dieu kien dung
        tmp4=0
        for i in range(0,vs_matrix.shape[0]):
            tmp4=tmp4+r1[i]**2
        tmp5=math.sqrt(tmp4)
        if tmp5 < 10**(-3):
            break
        #####
        beta= tmp4/tmp1   #beta
        p1= [r1[i]+beta*p0[i] for i in range(0,vs_matrix.shape[0])]
        p0=p1
        r0=r1
        x0=x1
        k=k+1
    elapsed_time_CGM = time.time() - start_time_CGM
    print("thoi gian giai he Ax=b:",elapsed_time_CGM)
    print(" so vong lap: ",k)
    return x1
def my_dblquad(two_var_function):
    result=0
    for i in range(0,7):
        result=result+two_var_function(GAUSS[i][0],GAUSS[i][1])*GAUSS[i][2]
    return result*0.5
def FEM_CGM(mesh,the_f_function):
    start_time_FEM_CGM = time.time()
    A = virtual_symmetric_matrix(mesh.count_iner_points())
    b = [0 for i in range(0,mesh.count_iner_points())]
    all_tri= mesh.all_triangle()
    for tri in all_tri:#J=(xj −xi)(yk −yi)−(xk −xi)(yj −yi)
        points=[mesh.points[tri[0]],mesh.points[tri[1]],mesh.points[tri[2]]]
        J=abs((points[1][0]-points[0][0])*(points[2][1]-points[0][1])-(points[2][0]-points[0][0])*(points[1][1]-points[0][1]))
        if tri[0] < A.shape[0]:
            A.set(tri[0],tri[0],A.get(tri[0],tri[0])+((points[1][1]-points[2][1])**2+(points[1][0]-points[2][0])**2)/(2*J))
            b[tri[0]]=b[tri[0]]+dblquad(lambda n,e: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*(1-e-n)*J, 0, 1, lambda e: 0, lambda e: 1-e)[0]
        if tri[1] < A.shape[0]:
            A.set(tri[1],tri[1],A.get(tri[1],tri[1])+((points[0][1]-points[2][1])**2+(points[0][0]-points[2][0])**2)/(2*J))
            b[tri[1]]=b[tri[1]]+dblquad(lambda n,e: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*e*J, 0, 1, lambda e: 0, lambda e: 1-e)[0]
        if tri[2] < A.shape[0]:
            A.set(tri[2],tri[2],A.get(tri[2],tri[2])+((points[1][1]-points[0][1])**2+(points[1][0]-points[0][0])**2)/(2*J))
            b[tri[2]]=b[tri[2]]+dblquad(lambda n,e: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*n*J, 0, 1, lambda e: 0, lambda e: 1-e)[0]
        if (tri[0] < A.shape[0]) and (tri[1] < A.shape[0]):
            A.set(tri[0],tri[1],A.get(tri[0],tri[1])+((points[1][1]-points[2][1])*(points[2][1]-points[0][1])+(points[1][0]-points[2][0])*(points[2][0]-points[0][0]))/(2*J))
        if (tri[1] < A.shape[0]) and (tri[2] < A.shape[0]):
            A.set(tri[1],tri[2],A.get(tri[1],tri[2])+((points[2][1]-points[0][1])*(points[0][1]-points[1][1])+(points[2][0]-points[0][0])*(points[0][0]-points[1][0]))/(2*J))
        if (tri[0] < A.shape[0]) and (tri[2] < A.shape[0]):
            A.set(tri[0],tri[2],A.get(tri[0],tri[2])+((points[0][1]-points[1][1])*(points[1][1]-points[2][1])+(points[0][0]-points[1][0])*(points[1][0]-points[2][0]))/(2*J))
    elapsed_time_FEM_CGM = time.time() - start_time_FEM_CGM
    print("thoi gian khoi tao A, b:",elapsed_time_FEM_CGM)
    return CGM(A, b)
def FEM_py(mesh,the_f_function):
    start_time_FEM_py = time.time()
    A = np.zeros((mesh.count_iner_points(),mesh.count_iner_points()))
    b = np.zeros(mesh.count_iner_points())
    all_tri= mesh.all_triangle()
    for tri in all_tri:#J=(xj −xi)(yk −yi)−(xk −xi)(yj −yi)
        points=[mesh.points[tri[0]],mesh.points[tri[1]],mesh.points[tri[2]]]
        J=abs((points[1][0]-points[0][0])*(points[2][1]-points[0][1])-(points[2][0]-points[0][0])*(points[1][1]-points[0][1]))
        if tri[0] < A.shape[0]:
            A[tri[0],tri[0]]=A[tri[0],tri[0]]+((points[1][1]-points[2][1])**2+(points[1][0]-points[2][0])**2)/(2*J)
            b[tri[0]]=b[tri[0]]+dblquad(lambda n,e: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*(1-e-n)*J, 0, 1, lambda e: 0, lambda e: 1-e)[0]
        if tri[1] < A.shape[0]:
            A[tri[1],tri[1]]=A[tri[1],tri[1]]+((points[0][1]-points[2][1])**2+(points[0][0]-points[2][0])**2)/(2*J)
            b[tri[1]]=b[tri[1]]+dblquad(lambda n,e: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*e*J, 0, 1, lambda e: 0, lambda e: 1-e)[0]
        if tri[2] < A.shape[0]:
            A[tri[2],tri[2]]=A[tri[2],tri[2]]+((points[1][1]-points[0][1])**2+(points[1][0]-points[0][0])**2)/(2*J)
            b[tri[2]]=b[tri[2]]+dblquad(lambda n,e: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*n*J, 0, 1, lambda e: 0, lambda e: 1-e)[0]
        if (tri[0] < A.shape[0]) and (tri[1] < A.shape[0]):
            A[tri[0],tri[1]]=A[tri[0],tri[1]]+((points[1][1]-points[2][1])*(points[2][1]-points[0][1])+(points[1][0]-points[2][0])*(points[2][0]-points[0][0]))/(2*J)
        if (tri[1] < A.shape[0]) and (tri[2] < A.shape[0]):
            A[tri[1],tri[2]]=A[tri[1],tri[2]]+((points[2][1]-points[0][1])*(points[0][1]-points[1][1])+(points[2][0]-points[0][0])*(points[0][0]-points[1][0]))/(2*J)
        if (tri[0] < A.shape[0]) and (tri[2] < A.shape[0]):
            A[tri[0],tri[2]]=A[tri[0],tri[2]]+((points[0][1]-points[1][1])*(points[1][1]-points[2][1])+(points[0][0]-points[1][0])*(points[1][0]-points[2][0]))/(2*J)
    for i in range(0,A.shape[0]):
        for j in range(0,i):
            A[i,j]=A[j,i]
    #print(A[:4],b)
    elapsed_time_FEM_py = time.time() - start_time_FEM_py
    print("thoi gian khoi tao A, b:",elapsed_time_FEM_py)
    start_time_py_sovle = time.time()
    resultfem = np.linalg.solve(A, b)
    elapsed_time_py_sovle = time.time() - start_time_py_sovle
    print("thoi gian giai he Ax=b:",elapsed_time_py_sovle)
    return resultfem
def FEM_CGM_GAUSS(mesh,the_f_function): # 'the_f_function' dc dinh nghia o tren
    start_time_FEM_CGM_GAUSS = time.time()
    A = virtual_symmetric_matrix(mesh.count_iner_points())
    b = [0 for i in range(0,mesh.count_iner_points())]
    all_tri= mesh.all_triangle()
    for tri in all_tri:#J=(xj −xi)(yk −yi)−(xk −xi)(yj −yi)
        points=[mesh.points[tri[0]],mesh.points[tri[1]],mesh.points[tri[2]]]
        J=abs((points[1][0]-points[0][0])*(points[2][1]-points[0][1])-(points[2][0]-points[0][0])*(points[1][1]-points[0][1]))
        if tri[0] < A.shape[0]:
            A.set(tri[0],tri[0],A.get(tri[0],tri[0])+((points[1][1]-points[2][1])**2+(points[1][0]-points[2][0])**2)/(2*J))
            b[tri[0]]=b[tri[0]]+my_dblquad(lambda e,n: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*(1-e-n)*J)
        if tri[1] < A.shape[0]:
            A.set(tri[1],tri[1],A.get(tri[1],tri[1])+((points[0][1]-points[2][1])**2+(points[0][0]-points[2][0])**2)/(2*J))
            b[tri[1]]=b[tri[1]]+my_dblquad(lambda e,n: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*e*J)
        if tri[2] < A.shape[0]:
            A.set(tri[2],tri[2],A.get(tri[2],tri[2])+((points[1][1]-points[0][1])**2+(points[1][0]-points[0][0])**2)/(2*J))
            b[tri[2]]=b[tri[2]]+my_dblquad(lambda e,n: the_f_function((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])*n*J)
        if (tri[0] < A.shape[0]) and (tri[1] < A.shape[0]):
            A.set(tri[0],tri[1],A.get(tri[0],tri[1])+((points[1][1]-points[2][1])*(points[2][1]-points[0][1])+(points[1][0]-points[2][0])*(points[2][0]-points[0][0]))/(2*J))
        if (tri[1] < A.shape[0]) and (tri[2] < A.shape[0]):
            A.set(tri[1],tri[2],A.get(tri[1],tri[2])+((points[2][1]-points[0][1])*(points[0][1]-points[1][1])+(points[2][0]-points[0][0])*(points[0][0]-points[1][0]))/(2*J))
        if (tri[0] < A.shape[0]) and (tri[2] < A.shape[0]):
            A.set(tri[0],tri[2],A.get(tri[0],tri[2])+((points[0][1]-points[1][1])*(points[1][1]-points[2][1])+(points[0][0]-points[1][0])*(points[1][0]-points[2][0]))/(2*J))
    elapsed_time_FEM_CGM_GAUSS = time.time() - start_time_FEM_CGM_GAUSS
    print("thoi gian khoi tao A, b:",elapsed_time_FEM_CGM_GAUSS)
    return CGM(A, b)
def draw_nghiem(mesh, nghiem):
    ax = a3.Axes3D(pl.figure())
    for tri in mesh.all_triangle():
        tmp= np.zeros((3,3))
        for i in range(0,3):
            for j in range(0,2):
                tmp[i,j]=mesh.points[tri[i]][j]
            if tri[i]<mesh.count_iner_points():
                tmp[i,2]=nghiem[tri[i]]     
        triagle = a3.art3d.Poly3DCollection([tmp])
        triagle.set_edgecolor('k')
        ax.add_collection3d(triagle)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    pl.show()
def draw_func(u,radius):
    R = np.linspace(0, radius, 1000)
    t = np.linspace(0,  2*np.pi, 1000)
    x = np.outer(R, np.cos(t))
    y = np.outer(R, np.sin(t))
    z = np.zeros(x.shape)
    for i in range(0,x.shape[0]):
        for j in range(0,x.shape[0]):
            z[i,j]= u(x[i,j],y[i,j])
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(x,y,z, cmap=cm.winter,
                           linewidth=0, antialiased=False)
    ax.set_zlim(-1,1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
def FEM_L2_error(mesh, c,nghiem_thuc):
##'c' la nghiem cua he FEM ung voi 'mesh', 'nghiem_thuc' dc dinh nghia o tren
    error=0
    c_size=mesh.count_iner_points()
    for tri in mesh.all_triangle():
        points=[mesh.points[tri[0]],mesh.points[tri[1]],mesh.points[tri[2]]]
        J=abs((points[1][0]-points[0][0])*(points[2][1]-points[0][1])-(points[2][0]-points[0][0])*(points[1][1]-points[0][1]))
        tmpi=tmpj=tmpk=0
        if tri[0]< c_size:
            tmpi= c[tri[0]]
        if tri[1]< c_size:
            tmpj= c[tri[1]]
        if tri[2]< c_size:
            tmpk= c[tri[2]]
        error=error + J*my_dblquad(lambda n,e: (nghiem_thuc((points[1][0]-points[0][0])*e+(points[2][0]-points[0][0])*n+points[0][0],(points[1][1]-points[0][1])*e+(points[2][1]-points[0][1])*n+points[0][1])-tmpi*(1-e-n)-tmpj*e-tmpk*n)**2)
    return math.sqrt(error)
def GiaTriNghiemTai(mesh,c,x,y):
    # xac dinh gia tri nghiem tai (x,y)
    # voi 'c' la nghiem cua FEM ung voi 'mesh'
    if math.sqrt(x**2+y**2)> mesh.radius:
        print("Error: (x,y) khong thuoc tap xac dinh")
    else:
        c_size=mesh.count_iner_points()
        tmpi=tmpj=tmpk=0
        for tri in mesh.all_triangle():
            points=[mesh.points[tri[0]],mesh.points[tri[1]],mesh.points[tri[2]]]
            s1=abs((points[1][0]-x)*(points[2][1]-y)-(points[2][0]-x)*(points[1][1]-y))
            s2=abs((x-points[0][0])*(points[2][1]-points[0][1])-(points[2][0]-points[0][0])*(y-points[0][1]))
            s3=abs((points[1][0]-points[0][0])*(y-points[0][1])-(x-points[0][0])*(points[1][1]-points[0][1]))
            J= abs((points[1][0]-points[0][0])*(points[2][1]-points[0][1])-(points[2][0]-points[0][0])*(points[1][1]-points[0][1]))
            if not(((s1+s2) >J) or ((s1+s3) >J) or ((s3+s2) >J)):
                e=((points[2][1]-points[0][1])*(x-points[0][0])-(points[2][0]-points[0][0])*(y-points[0][1]))/J
                n=(-(points[1][1]-points[0][1])*(x-points[0][0])+(points[1][0]-points[0][0])*(y-points[0][1]))/J
                if tri[0]< c_size:
                   tmp_i= c[tri[0]]
                if tri[1]< c_size:
                   tmp_j= c[tri[1]]
                if tri[2]< c_size:
                   tmp_k= c[tri[2]]
                return tmp_i*(1-e-n) + tmp_j*e +tmp_k*n
    print("Error: (x,y) khong thuoc tap xac dinh")
                
            
start_time_total= time.time()
start_time_mesh = time.time()
M=cmesh.cirle_mesh(_resolution,_Radius)
elapsed_time_mesh  = time.time() - start_time_mesh 
print("thoi gian khoi tao mesh:",elapsed_time_mesh)
print("    + Do phan giai:",_resolution)
print("    + So diem trong cua mesh (khong tinh dinh tai bien):",M.count_iner_points())
print("    + So phan tu huu han ( so tam giac phan hoach ):",M.num_Rec())
c =FEM_CGM_GAUSS(M,_f)
elapsed_time_total = time.time() - start_time_total
print("Tong thoi gian:",elapsed_time_total )
print("h (h = R/resolution) : ", M.radius/M.size)
print("Sai so cua nghiem tinh theo chuan L2 :", FEM_L2_error(M, c, _u))
u_x_y=GiaTriNghiemTai(M,c,_toa_do[0],_toa_do[1])
print("Gia tri cua nghiem tai ",_toa_do,"la:",u_x_y, "      Error:",abs(_u(_toa_do[0],_toa_do[1])-u_x_y))
if _draw_nghiem:
    draw_nghiem(M,c)
if _draw_u:
    draw_func(_u,M.radius)
