# SPDX-License-Identifier: GPL-3.0-only
#
set(ID_DOC_INSTALL_DIR doc)

function(id_install_doc)
    cmake_parse_arguments(id_install "" "" "FILES" ${ARGV})
    install(FILES ${id_install_FILES} DESTINATION ${ID_DOC_INSTALL_DIR})
endfunction()

id_install_doc(FILES known-bugs.md debug-flags.md)
