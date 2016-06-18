TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXXFLAGS_RELEASE += -DARMA_NO_DEBUG -O3

SOURCES += main.cpp

HEADERS += kmcrngzig.h \
           zigrandom.h \
           zignor.h

SOURCES += zigrandom.cpp \
           zignor.cpp


include(deployment.pri)
qtcAddDeployment()

