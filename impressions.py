import numpy as np
import matplotlib.pyplot as plt

# учет миграционных процессов в первой и второй группе
V1 = 0.5
V2 = 0.3

# группы художников
u10 = 0.8  # следующие классическим канонам
u20 = 0.4  # новое течение

a = 0.1  # нижняя граница численности первой группы
B = 0.75  # при повышении u2 снижается
b = 0.75

m1 = 0.1  # удельная скорость

sigma11 = 3
sigma12 = 1
sigma21 = 1
sigma22 = 6

l = 10
N = 10  # кол-во отрезков разбиения отрезка [0, l]
h = l / N
#print("h = ",h)

tmax = 80
M = 10  # кол-во отрезков разбиения [0, tmax]
tao = tmax / M
#print("tao = ", tao)

x = np.linspace(0., l, N)
t = np.linspace(0., tmax, M)
#print("x = ", x)
#print("t = ", t)

u1 = [[0] * N for i in range(M)]
u2 = [[0] * N for i in range(M)]
for i in range(N):
    u1[0][i] = u10
for i in range(N):
    u2[0][i] = u20

# Вывод матрицы на экран
def print_arr( string, namevec, a ):
    if (type(a) == int) or (type(a) == float):
        print(a)
    else:
        print( string )
        for k in range(len(a)):
            print("{}[{}] = {:8.4f}".format(namevec, k, a[k]))


# Проверка 3х-диаг. матрицы коэффициентов на корректность
def isCorrectArray(a):
    n = len(a)

    for row in range(0, n):
        if( len(a[row]) != n ):
            print('Не соответствует размерность')
            return False

    for row in range(1, n - 1):
        if(abs(a[row][row]) < abs(a[row][row - 1]) + abs(a[row][row + 1])):
            print('Не выполнены условия достаточности')
            return False

    if (abs(a[0][0]) < abs(a[0][1]))or(abs(a[n - 1][n - 1]) < abs(a[n - 1][n - 2])):
        print('Не выполнены условия достаточности')
        return False


    for row in range(0, len(a)):
        if( a[row][row] == 0 ):
            print('Нулевые элементы на главной диагонали')
            return False
    return True


# Процедура нахождения решения 3-х диагональной матрицы
def solution(a, b):
    if( not isCorrectArray(a) ):
        print('Ошибка в исходных данных')
        return -1

    n = len(a)
    x = [0 for k in range(0, n)] # обнуление вектора решений
    print('Размерность матрицы: ',n,'x',n)

    # Прямой ход
    p = [0 for k in range(0, n)]
    q = [0 for k in range(0, n)]
    # для первой 0-й строки
    p[0] = a[0][1] / (-a[0][0])
    q[0] = ( - b[0]) / (-a[0][0])
    for i in range(1, n - 1): # заполняем за исключением 1-й и (n-1)-й строк матрицы
        p[i] = a[i][i+1] / ( -a[i][i] - a[i][i-1]*p[i-1] )
        q[i] = ( a[i][i-1]*q[i-1] - b[i] ) / ( -a[i][i] - a[i][i-1]*p[i-1] )
    # для последней (n-1)-й строки
    p[n-1] = 0
    q[n-1] = (a[n-1][n-2]*q[n-2] - b[n-1]) / (-a[n-1][n-1] - a[n-1][n-2]*p[n-2])

    print_arr('Прогоночные коэффициенты p: ','p', p)
    print_arr('Прогоночные коэффициенты q: ','q', q)

    # Обратный ход
    x[n-1] = q[n-1]
    for i in range(n-1, 0, -1):
        x[i-1] = p[i-1] * x[i] + q[i-1]

    return x

print("k=0")
for k in range(M):
  print(u1[k])

for k in range(M):
  print(u2[k])

for k in range(1, M):
  print("k=",k)

  # считаем коэф-ты
  A1 = np.zeros(N)
  B1 = np.zeros(N)
  C1 = np.zeros(N)
  F1 = np.zeros(N)

  A1[0] = V1/(h*h)
  B1[0] = -(2*V1)/(h*h) + m1*(u1[k-1][0] - a)*(1 - u1[k-1][0]) - B*u2[k-1][0] - 1/tao
  C1[0] =  V1/(h*h)
  F1[0] = -sigma11/tao

  A2 = np.zeros(N)
  B2 = np.zeros(N)
  C2 = np.zeros(N)
  F2 = np.zeros(N)

  A2[0] = V2/h*h
  B2[0] = -(2 * V2) / (h*h) - B * (b - u1[k-1][0]) - 1/tao
  C2[0] = V2/(h*h)
  F2[0] = -sigma21/tao

  for j in range(1, N-1):
    A1[j] = 0
    C1[j] = V1/(h*h)
    B1[j] = -(2*V1)/(h*h) + m1*(u1[k-1][j] - a)*(1 - u1[k-1][j]) - B*u2[k-1][j] - 1/tao
    F1[j] = -u1[k-1][j]/tao - V1*sigma11/h*h


    A2[j] = 0
    B2[j] = - 1/tao - B*(b - u1[k-1][j]) - 2*V2/h*h
    C2[j] = V2/(h*h)
    F2[j] = -u2[k-1][j]/tao - (sigma21*V2)/h*h

  A1[N-1] = V1/(h*h)
  B1[N-1] =  -(2*V1)/(h*h) + m1*(u1[k-1][j] - a)*(1 - u1[k-1][j]) - B*u2[k-1][j] - 1/tao
  C1[N-1] = 0
  F1[N-1] = -u1[k-1][j]/tao - (V1*sigma12)/h*h


  A2[N-1] = V2/(h*h)
  B2[N-1] = - 1/tao - B*(b - u1[k-1][j]) - 2*V2/h*h
  C2[N-1] = 0
  F2[N-1] = -u2[k-1][j]/tao - (sigma22*V2)/h*h

  #строим матрицу и считаем u1[k][j] j = 1...N
  Matrix1  =  [[0] * N for i in range(N)]
  Matrix2  =  [[0] * N for i in range(N)]

  for i in range(N):
    Matrix1[i][i] = B1[i]
    if(A1[i] != 0): Matrix1[i][i-1] = A1[i]
    if(C1[i] != 0): Matrix1[i][i+1] = C1[i]

  for i in range(N):
    Matrix2[i][i] = B2[i]
    if(A2[i] != 0): Matrix2[i][i-1] = A2[i]
    if(C2[i] != 0): Matrix2[i][i+1] = C2[i]


  for k in range(1, M):
    print("k=",k)

# ГРАФИКИ
  print("A1 = ",A1)
  print("B1 = ",B1)
  print("C1 = ",C1)
  print("F1 = ", F1)
  for i in range(N):
      print("Matrix1 = ", Matrix1[i])

  print("A2 = ",A2)
  print("B2 = ",B2)
  print("C2 = ",C2)
  print("F2 = ", F2)
  for i in range(N):
      print("Matrix2 = ", Matrix2[i])

  # метод прогонки
  print(isCorrectArray(Matrix1))
  print(isCorrectArray(Matrix2))

  x1 = solution(Matrix1, F1)  # Вызываем процедуру решения
  print_arr('Решение: ','x1', x1)

  x2 = solution(Matrix2, F2)  # Вызываем процедуру решения
  print_arr('Решение: ','x2', x2)

  for i in range(N):
    u1[k][i] = x1[i]


  for i in range(N):
    u2[k][i] = x2[i]


  print("u1[K] = ",u1[k])
  print("u2[K] = ",u2[k])

  plt.plot (t, u1[k], label = "u1")
  plt.plot (t, u2[k], label = "u2")
  #plt.xlabel('t')
  plt.grid(True)
  plt.legend()
  plt.show()

  fig1 = plt.figure()
  xtu1 = plt.axes(projection='3d')
  xtu1.plot3D(x, t, u1[k])
  xtu1.plot3D(x, t, u2[k])
  #plt.xlabel('x')
  #plt.ylabel('t')

  fig2 = plt.figure()
  xtu2 = plt.axes(projection='3d')
  xtu2.plot3D(x, t, u2[k], 'red')
  plt.xlabel('x')
  plt.ylabel('t')


  plt.show()