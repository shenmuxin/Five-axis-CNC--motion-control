import numpy as np
import sympy

DATA_PATH = './path_interpolation_xyzijk.txt'
SAVE_PATH = './path_interpolation_xyzbc.txt'

def read_path_data(file_path):
    """
    :@func: 数据读取函数
    :param file_path: 数据存放路径
    :return : 初始插值点
              shape(n, 6)
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    path_data = []
    for line in lines[1:]:
        line_data = [float(x) for x in line.strip().split(',')]
        path_data.append(line_data)
    return np.array(path_data)


def Rodrigues(rotation_vector, rotation_theta):
    """
    :@func: 计算旋转矩阵
    :param rotation_vector: 旋转向量
    :return : 
    """
    k = rotation_vector
    k_ = sympy.Matrix([[0, -k[2, 0], k[1, 0]], [k[2, 0], 0, -k[0, 0]], [-k[1, 0], k[0, 0], 0]])
    Rodrigues_mat = sympy.cos(rotation_theta) * sympy.eye(3) + (1 - sympy.cos(rotation_theta)) * k * k.T + sympy.sin(
        rotation_theta) * k_
    return Rodrigues_mat

def inv_kinema(path_data):
    """
    :@func: 求解逆运动学
    :param path_data: 刀心路径点信息 shape(num,6) (x,y,z,i,j,k)
    :return :  求解的刀轴信息 shape(num,5) (x,y,z,b,c)
    """
    #加载数据
    points = read_path_data(DATA_PATH)
    curve_points = points[:,0:3]    # 路径
    normal_vectors = points[:,3:]  # 法向量


    # 旋转轴初始方向
    wc = sympy.Matrix([[0], [0], [1]])
    wb = sympy.Matrix([[0], [0.5 * sympy.sqrt(2)], [0.5 * sympy.sqrt(2)]])

    theta_b, theta_c ,x,y,z= sympy.symbols('theta_b theta_c x y z')
    # 对C旋转
    cRodrigues = Rodrigues(wc, theta_c)
    # 对B旋转
    bRodrigues = Rodrigues(wb, theta_b)
    # print(bRodrigues)
    zero=np.zeros((3,1))
    # 生成旋转矩阵
    ec=np.append(cRodrigues,zero,axis=1)
    e_bu=np.array([0,0,0,1]).reshape(1,4)
    e_c=np.append(ec,e_bu,axis=0).reshape(4,4)
    eb=np.append(bRodrigues,zero,axis=1)
    e_b=np.append(eb,e_bu,axis=0).reshape(4,4)
    # print(eb,e_b)
    e_x=np.array([1,0,0,x,0,1,0,0,0,0,1,0,0,0,0,1]).reshape(4,4)
    e_y=np.array([1,0,0,0,0,1,0,y,0,0,1,0,0,0,0,1]).reshape(4,4)
    e_z=np.array([1,0,0,0,0,1,0,0,0,0,1,z,0,0,0,1]).reshape(4,4)
    # 初始位型
    mt0 = np.array([1,0,0,10,0,1,0,20,0,0,1,100,0,0,0,1]).reshape(4,4)
    gmw=np.eye(4)
    gmt1=np.dot(e_x,e_y)
    gmt2= np.dot(gmt1,e_z)
    gmt3=np.dot(e_c,e_b)
    gmt4=np.dot(gmt3,mt0)
    gmt=np.dot(gmt2,gmt4)
    # 转移矩阵
    gm=np.dot(gmw,gmt)

    # 带入解析式，求解
    bx = normal_vectors[:,0]/np.linalg.norm(normal_vectors[:,0])#归一化
    by = normal_vectors[:,1]/np.linalg.norm(normal_vectors[:,1])
    bz = normal_vectors[:,2]/np.linalg.norm(normal_vectors[:,2])
    x = curve_points[:,0]
    y = curve_points[:,1]
    z = curve_points[:,2] 

    print(2*bz-1)
    theta_b = np.arccos(2*bz-1)
    theta_c = np.arcsin(((bz-1)*bx+np.sqrt(2)*np.sqrt(bz-bz**2)*by)/(1-bz**2))

    # 存数据
    myfile = open("path_interpolation_xyzbc.txt", "w")
    myfile.write('x   y   z   b   c\n')
    for i in range(0, normal_vectors.__len__()):
        myfile.write('{:>6.3f} '.format(x[i]))
        myfile.write('{:>6.3f} '.format(y[i]))
        myfile.write('{:>6.3f} '.format(z[i]))
        myfile.write('{:>6.3f} '.format(theta_b[i]))
        myfile.write('{:>6.3f}\n'.format(theta_c[i]))

    myfile.close()