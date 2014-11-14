QT += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
CONFIG -= app_bundle
CONFIG += opengl
TEMPLATE = app
CONFIG += console
CONFIG(debug, debug | release) {
  MODE = "debug"
} else {
  MODE = "release"
}

win32|win64 {
  #TMP_FILE_PREFIX = "D:/tmp/$$PWD"
  TMP_FILE_PREFIX = "C:/tmp/$$MODE"
  DEFINES += 'CURRENT_DIRECTORY=$$PWD/data/'
  QMAKE_CXXFLAGS += /O2
} else {
  TMP_FILE_PREFIX = "/tmp/$$PWD/$$MODE"
  QMAKE_CXXFLAGS += -fopenmp
  QMAKE_CXXFLAGS_RELEASE -= -O2
  QMAKE_CXXFLAGS_RELEASE += -O3
  QMAKE_LFLAGS_RELEASE += -O3
  QMAKE_LFLAGS_RELEASE -= -O1
  QMAKE_CXXFLAGS_DEBUG += -g
  DEFINES += 'CURRENT_DIRECTORY=\'\"$$PWD/data/\"\''
  DEFINES += 'DATA_DIRECTORY=\'\"$$PWD/data/\"\''
}

macx {
  CC_VERSION = 4.8
  QMAKE_CC = /opt/local/bin/gcc-mp-$$CC_VERSION
  QMAKE_CXX = /opt/local/bin/g++-mp-$$CC_VERSION
  QMAKE_LINK = /opt/local/bin/g++-mp-$$CC_VERSION
  QMAKE_CXXFLAGS += -std=c++11
  QMAKE_LFLAGS += -fopenmp
  #QMAKE_CXX = clang++
  LIBS+= -framework GLUT -framework OpenGL -L/opt/local/lib -lglew
  #LIBS += -lgslcblas -larmadillo
  INCLUDEPATH += /opt/local/include
  INCLUDEPATH += /usr/include
} else:unix {
  QMAKE_CC = gcc-4.8
  QMAKE_CXX = g++-4.8
  QMAKE_LINK = g++-4.8
  QMAKE_CXXFLAGS += -std=c++11  -Wno-unused-result
  LIBS += -lglut -lGL -lGLU -lGLEW
  #LIBS += -lpthread -fopenmp
  #LIBS += -lgslcblas  -larmadillo
} else:win32|win64 {
  LIBS += -L"C:/lib/lib/x64"
  LIBS += GlU32.lib freeglut.lib glew32.lib OpenGL32.lib
  #LIBS += libgslcblas.a
  INCLUDEPATH += "C:/lib/include"
}

msvc {
  #msvc:QMAKE_CXXFLAGS += -openmp
  QMAKE_CXXFLAGS += /openmp
  DEFINES += _CRT_SECURE_NO_WARNINGS
}

OTHER_FILES += \
    $$PWD/data/scene.conf \
    $$PWD/data/global.conf \
    ../Dropbox/pbd-rod/CMakeLists.txt

#TARGET = skeleton
#DESTDIR = $$PWD
#OBJECTS_DIR = $$TMP_FILE_PREFIX/$$TARGET/obj
#MOC_DIR = $$TMP_FILE_PREFIX/$$TARGET/moc
#RCC_DIR = $$TMP_FILE_PREFIX/$$TARGET/rcc
#UI_DIR = $$TMP_FILE_PREFIX/$$TARGET/gui

include($$PWD/src/lib/lib.pri)
include($$PWD/src/rod.pri)
