INCLUDEPATH += $$PWD

SOURCES += \
    $$PWD/main.cpp\
    $$PWD/main_window.cpp \
    $$PWD/open_gl_qt.cpp \
    $$PWD/global.cpp \
    $$PWD/pbd_elastic_rod.cpp \
    $$PWD/simulator.cpp \
    $$PWD/elastic_rod_shape.cpp \
    $$PWD/frame_attachment.cpp \

HEADERS  += \
    $$PWD/main_window.h \
    $$PWD/open_gl_qt.h \
    $$PWD/pbd_elastic_rod.h \
    $$PWD/simulator.h \
    $$PWD/elastic_rod_shape.h \
    $$PWD/frame_attachment.h \
    $$PWD/global.h \

FORMS += \
    $$PWD/ui/main_window.ui
