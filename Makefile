# (1)コンパイラ
CC  = g++
# (2)コンパイルオプション
# CFLAGS    =
# (3)実行ファイル名
TARGET  = Pso
# (4)コンパイル対象のソースコード
SRCS    += main.cpp \
		   pso.cpp

#HEADERS += pso.h

# (5)オブジェクトファイル名
OBJS    = $(SRCS:.cpp=.o)
 
# (6)インクルードファイルのあるディレクトリパス
#INCDIR += 

# (7)ライブラリファイルのあるディレクトリパス

LIBDIR += -I/usr/local/Cellar/eigen/3.3.8_1/include/
#LIBDIR += /Users/maru/Downloads/boost_1_76_0
#LIBDIR += /usr/local/Cellar/boost/1.75.0_2/include/boost/

# (8)追加するライブラリファイル
LIBS += -L/usr/local/Cellar/eigen/3.3.8_1/include/eigen3/Eigen
#LIBS += -L/Users/maru/Downloads/boost_1_76_0/stage/lib -lboost_system -lboost_filesystem
#LIBS += -L/Users/maru/Downloads/boost_1_76_0 -lboost_system -lboost_filesystem
#LIBS += -L/usr/local/Cellar/boost/1.75.0_2/include/boost -lboost_system -lboost_filesystem
#LIBS += -L/usr/local/Cellar/boost/1.75.0_2/include/boost/
# (9)ターゲットファイル生成
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBDIR) $(LIBS)
	
# (10)オブジェクトファイル生成
$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(INCDIR) -c $(SRCS)

# (11)"make all"で make cleanとmakeを同時に実施。
all: clean $(OBJS) $(TARGET)
# (12).oファイル、実行ファイル、.dファイルを削除
clean:
	-rm -f $(OBJS) $(TARGET) *.d