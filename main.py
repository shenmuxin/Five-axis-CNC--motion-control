import numpy as np
import matplotlib.pyplot as plt
from kinematics_inverse import inv_kinema, read_path_data

#定义全局常量
MAX_SPEED = 0.05           #最大合成速度 m/s
MAX_ACCEL = 0.5             #最大合成加速度 m/2^2
CORNER_TIME = 0.003         #拐弯时间 s
PERIOD = 0.001              #插补周期 s,


def write_file(file_path, data, string=''):
    """
    :@func: 数据保存函数
    :param file_path: 数据存放路径
    :return : 
    """
    with open(file_path, 'w') as f:
        f.write(string + '\n')
        for i in range(len(data)):
            f.write(','.join(list(map(str, data[i]))) + '\n')

def write2file(file_path, data, string=''):
    with open(file_path, 'w') as f:
        f.write(string + '\n')
        for i in range(len(data)):
            f.write('{:>6.3f}\n'.format(data[i]))


def plot_trajectory(traj_task, origin_traj):
    """
    :@func: 绘制仿真图形
    :param traj_task:   小线段速度规划得到的插补点
    :param origin_traj: 初始控制点
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.plot([p[0] for p in traj_task], [p[1] for p in traj_task], [p[2] for p in traj_task], 'b')
    ax.plot([p[0] for p in traj_task], [p[1] for p in traj_task], [p[2] for p in traj_task], c='r')
    ax.scatter([p[0] for p in origin_traj], [p[1] for p in origin_traj], [p[2] for p in origin_traj], 'k')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


def speed_planning(traj_data, max_speed=MAX_SPEED, max_accel=MAX_ACCEL, corner_time=CORNER_TIME, period=PERIOD):
    """
    :@func:  参考《An optimal feedrate model and solution algorithm for a 
                    high-speed machine of small line blocks with look-ahead》给出的方法,进行速度规划
    :param traj_data: 路径信息shape(num,3) (x,y,z)
    :param max_speed=0.05: 最大合成速率 m/s
    :param max_accel=0.5: 最大合成加速度 m/s^2
    :param corner_time=0.003: 拐弯时间 s
    :param period = 0.001: 插补周期 s
    :return:
    """
    # 定义常量
    block_max = 4   # 最大前瞻数
    num = len(traj_data)
    # 开辟存储空间
    a = np.zeros([num])
    b = np.zeros([num])
    e = np.zeros([num])
    c = np.zeros([num])
    d = np.power(max_speed,2)
    Vs = np.zeros([num,1])
    # 节点差分
    dp = np.diff(traj_data, axis=0)
    # 初始化
    e[0] = d
    b[0] = 0
    b[1] = 2*max_accel*np.sqrt(np.power(dp[0],2).sum())

    for i in range(1,num-1):  # Vs中第一个速度和最后一个速度都为0
        j = 0
        if j < block_max:
            if i+j <num:
                j=j+1
                b[i+j] = 2*max_accel*np.sqrt(np.power(dp,2).sum())
                e[i+j-1] = corner_time*max_accel*max_accel / (2*(1-(dp[i+j-2,0]*dp[i+j-1,0]
                                                                    +dp[i+j-2,1]*dp[i+j-1,1]
                                                                    +dp[i+j-2,2]*dp[i+j-1,2])
                                                                    /(np.sqrt(np.power(dp[i+j-2],2).sum())
                                                                      *np.sqrt(np.power(dp[i+j-1],2).sum()))))
                c[i+j-1] = np.min([e[i+j-1],d])
                if c[i+j-1] > b[i+j]:
                    sum = 0
                    for k in range(i+1,i+j):
                        sum += b[k]
                    if sum >= c[i]:
                        Vs[i]=np.sqrt(np.min([Vs[i-1]*Vs[i-1]+b[i], e[i]]))
                elif c[i+j-1] <= b[i+j]:
                    sum = 0
                    for k in range(i+1,i+j):
                        sum += b[k]
                    Vs[i] = np.sqrt(np.min([sum+c[i+j-1], Vs[i-1]*Vs[i-1]+b[i], e[i]]))

            elif i+j >= num:
                sum = 0
                for k in range(i+1,i+j):
                    sum += b[k]
                Vs[i] = np.sqrt(np.min([sum, Vs[i-1]*Vs[i-1]+b[i], e[i]]))

        elif j>= block_max:
            sum = 0
            for k in range(i+1,i+j):
                sum += b[k]
                Vs[i] = np.sqrt(np.min([sum, Vs[i-1]*Vs[i-1]+b[i], e[i]]))
    return Vs


def calc_lines_info(path_data, plan_vels, max_speed=MAX_SPEED, max_accel=MAX_ACCEL):
    """
    :@func: 三轴小线段计算中间信息函数
    :param path_data:  初始路径信息 shape(n,6) (x,y,z,i,j,k)
    :param plan_vels:  初始路径点的合成规划速度 shape(n,1)
    :param max_speed = 0.05: 三轴最大合成速度  m/s
    :param max_accel = 0.5: 三轴最大合成加速度 m/s^2
    :return : 

    """
    
    num = path_data.shape[0]
    # distances = np.zeros([num -1, 1])  #插值点间的距离shape(n-1, 3)
    dp = np.diff(path_data[:,0:3],axis=0)
    # 开辟存储空间
    vels_m = np.zeros([num -1])
    s1 = np.zeros([num -1])
    s2 = np.zeros([num -1])
    s3 = np.zeros([num -1])
    ta = np.zeros([num -1])
    td = np.zeros([num -1])
    tl = np.zeros([num -1])
    for i in range(num -1):
        # 计算距离
        # distances[i] = np.sqrt(np.power(path_data[i+1,0:3] - path_data[i,0:3],2).sum()) 
        dis_i = np.sqrt(np.power(dp[i],2).sum())
        # 计算合成速度
        vi = plan_vels[i]
        vi_1 = plan_vels[i+1]
        
        # vels_m[i] = min(np.sqrt((np.power(vi,2) + np.power(vi_1,2) + 2*max_accel*distances[i])/2), max_speed)
        vels_m[i] = min(np.sqrt((np.power(vi,2) + np.power(vi_1,2) + 2*max_accel*dis_i)/2), max_speed)
        s1[i] = (np.power(vels_m[i],2) - np.power(vi,2))/2
        s3[i] = (np.power(vels_m[i],2) - np.power(vi_1,2))/2
        s2[i] = dis_i - s1[i] -s3[i]
        ta[i] = (vels_m[i] - vi) / max_accel
        td[i] = (vels_m[i] - vi_1) / max_accel
        tl[i] = s2[i] / vels_m[i]
    
    return vels_m, ta, td, tl


def calc_axis_point(start, end, Vi, Vi_1, Vm, timea, timed, timel, max_accel=MAX_ACCEL, period=PERIOD):
    """
    :@func:  计算一段的插补点
    :param start:  shape(6,)
    :param end: shape(6,)
    :param Vm: 该段的最大速度
    :param Vi: 起始速度
    :param Vi_1: 终点速度
    :param timea: 加速时间
    :param timed: 匀速时间
    :param timel: 减速时间

    :return: 返回给定两点之间的插补点,shape(n,6)
    """

    dist = np.sqrt(np.power(end[0:3]-start[0:3],2).sum())
    # 计算三轴x,y,z分别所占的比例
    k1 = (end[0:3]-start[0:3])[0] / np.sqrt(np.power(end[0:3]-start[0:3],2).sum())
    k2 = (end[0:3]-start[0:3])[1] / np.sqrt(np.power(end[0:3]-start[0:3],2).sum())
    k3 = (end[0:3]-start[0:3])[2] / np.sqrt(np.power(end[0:3]-start[0:3],2).sum())

    # 计算法向量的比例
    delta = (end[3:]-start[3:]) / (timea+timed+timel)

    # 插值周期数
    t_count = 1
    V1 = Vi         #起始点速度
    V2 = Vi_1       #终点速度

    # 开辟存储空间
    line_points = np.zeros([1,6])
    dis = 0
    while(t_count*period <= timea+timed+timel and dis<=dist):
        # 计算加速时间内的插值点
        if t_count*period <= timea:
            dis += V1*period + 0.5*max_accel*period*period
            V1 += max_accel*period
            if V1 > Vm:
                V1 = Vm
            t_count += 1
            delta_vec = delta * t_count * period
            delta_xyz = np.concatenate([k1*dis, k2*dis, k3*dis], axis=0)
            dp = np.concatenate([delta_xyz, delta_vec], axis=0)
            p = start + dp
            line_points = np.vstack([line_points, p])

        # 计算匀速时间内的插值点
        elif timea < t_count*period <= timea+timel and timel!=0:
            dis += V1 * period
            t_count += 1
            delta_vec = delta * t_count * period
            delta_xyz = np.concatenate([k1*dis, k2*dis, k3*dis], axis=0)
            dp = np.concatenate([delta_xyz, delta_vec], axis=0)
            p = start + dp
            line_points = np.vstack([line_points, p])
        
        # 计算减速时间内的插值点
        elif timea+timel < t_count*period <= timea+timel+timed:
            dis += V1*period - 0.5*max_accel*period*period
            V1 = V1 - max_accel*period
            if V1 < V2:
                V1 = V2
            t_count += 1
            delta_vec = delta * t_count * period
            delta_xyz = np.concatenate([k1*dis, k2*dis, k3*dis], axis=0)
            dp = np.concatenate([delta_xyz, delta_vec], axis=0)
            p = start + dp
            line_points = np.vstack([line_points, p])
    
    return line_points[1:]



def calc_axis_points(path_data, plan_vels, max_speed=MAX_SPEED, max_accel=MAX_ACCEL, corner_time=CORNER_TIME, period=PERIOD):
    """
    :@func:  计算所有的插补点
    :param path_data:  shape(num,6)
    :param plan_vels: shape(num,)
    :param max_speed = 0.05: 三轴最大合成速度  m/s
    :param max_accel = 0.5: 三轴最大合成加速度 m/s^2
    :param corner_time = 0.003: 拐弯时间 s
    :param period = 0.001: 固定插补周期 s

    :return: 所有的插补点,shape(n,6)
    """
    vels_m, ta, td, tl = calc_lines_info(path_data, plan_vels, max_speed, max_accel)

    # 开辟存储空间
    axis_points = np.zeros([1,6])
    num = len(path_data) - 1
    for i in range(num):
        line_points = calc_axis_point(path_data[i], path_data[i+1], plan_vels[i], plan_vels[i+1], vels_m[i], ta[i], td[i], tl[i], max_accel, period)
        axis_points = np.vstack([axis_points, line_points])
    
    return axis_points[1:]


def diff_vel_accel(axis_points, period=PERIOD):
    """
    :@func:  通过五轴插补点计算实际速度和实际加速度
    :param axis_points:  实际插补点
    :param period: 插补周期

    :return: 实际速度, 实际加速度
    """

    delta_d = np.diff(axis_points[:,0:3], axis=0)
    dv = delta_d / period
    real_vels = np.sqrt(np.power(dv,2).sum(axis=1))
    delta_v = np.diff(dv, axis=0)
    da = delta_v / period
    real_accels = np.sqrt(np.power(da,2).sum(axis=1))
    # print(real_vels[0:5])
    # print(real_accels[0:5])

    # 将实际速度写入到文件real_velocities.txt中
    write2file('real_velocities.txt', real_vels, string = "real velocities")
    # 将实际加速度写入到文件real_accelerations.txt中
    write2file('real_accelerations.txt', real_accels, string="real accelerations")
    return real_vels, real_accels


def sliding_average(data, window_size):
    """
    :@func: 实现滑动平均滤波
    :param data: 滑动滤波的数据
    :param window_size:  滑动窗口大小
    """
    filtered_data = []
    for i in range(len(data)):
        if i < window_size:
            filtered_data.append(sum(data[:i+1]) / (i+1))
        else:
            filtered_data.append(sum(data[i-window_size+1:i+1]) / window_size)
    
    write2file("filtered_accelerations.txt", filtered_data, string = "filtered accelerations")
    return filtered_data

    


if __name__ == "__main__":

    # 读取HW3所得到的路径点
    path_data = read_path_data('original_path.txt')

    # 进行速度规划,得到初始路径点上的规划速度
    plan_vels = speed_planning(path_data[:,0:3])
    # 保存规划速度到planned_velocities.txt中
    write_file('planned_velocities.txt', plan_vels,'')

    axis_points = calc_axis_points(path_data, plan_vels)
    print('\n插补点的个数为:', len(axis_points),'\n')

    # 保存插补点到path_interpolation.txt中
    write_file('path_interpolation_xyzijk.txt', axis_points, string='x   y   z   i   j   k')

    # 逆运动学求解
    inv_kinema(axis_points)

    # 求解实际的速度和加速度,保存到real_velocities.txt和real_accelerations.txt中
    _, real_accels = diff_vel_accel(axis_points)

    # 对实际加速度进行滑动滤波
    filtered_accels = sliding_average(real_accels, 4)

    # 绘制生成的路径图
    plot_trajectory(axis_points, path_data)











