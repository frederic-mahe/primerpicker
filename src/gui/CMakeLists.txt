# The graphical front-end (Qt 6 Widgets): a library (also linked by the
# GUI tests) and the program. Presentation only: every value is
# computed by the core library.
add_library(primerpicker_gui
  analysis_view.cpp
  direct_entry_dialog.cpp
  main_window.cpp)
set_target_properties(primerpicker_gui PROPERTIES AUTOMOC ON)
target_link_libraries(primerpicker_gui PUBLIC primerpicker_core Qt6::Widgets)
target_compile_options(primerpicker_gui PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
# QString conversions and Qt's keyword macros are explicit; no API
# deprecated in the minimum Qt version
target_compile_definitions(primerpicker_gui PUBLIC
  QT_NO_CAST_FROM_ASCII
  QT_NO_CAST_TO_ASCII
  QT_NO_CAST_FROM_BYTEARRAY
  QT_NO_KEYWORDS
  QT_NO_NARROWING_CONVERSIONS_IN_CONNECT
  QT_NO_URL_CAST_FROM_STRING
  QT_DISABLE_DEPRECATED_UP_TO=0x060800)

add_executable(primerpicker-gui WIN32 MACOSX_BUNDLE
  main.cpp)
target_link_libraries(primerpicker-gui PRIVATE primerpicker_gui)
target_compile_options(primerpicker-gui PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
