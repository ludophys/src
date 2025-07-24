project( 'neutron_bg'
       , 'cpp'
       , version : 'v0.0.0'
       , default_options : [ 'buildtype=debugoptimized'
                           , 'debug=true'
                           , 'optimization=2'
                           , 'cpp_std=c++20'
                           , 'prefix=@0@/../install/neutron_bg/'.format(meson.source_root())
                           , 'libdir=lib'
                           ]
       )

nain4  = dependency( 'nain4'
                   , method  : 'pkg-config'
                   , required: true
                   )

geant4 = dependency( 'Geant4'
                   , method     : 'cmake'
                   , required   : true
                   , components : ['ui_all', 'vis_all']
                   , modules    : nain4.get_variable(pkgconfig: 'geant4_modules').split(','))

neutron_bg_lib_deps     = [nain4, geant4]
neutron_bg_lib_include  = include_directories('.')
neutron_bg_lib_includes = ['LXe.hh']
neutron_bg_lib_sources  = ['LXe.cc']

geant4_include = geant4.get_variable(cmake    : 'Geant4_INCLUDE_DIRS')
nain4_include  =  nain4.get_variable(pkgconfig: 'includedir'         )


neutron_bg_lib = shared_library( 'neutron_bg'
                                   , neutron_bg_lib_sources
                                   , include_directories: [neutron_bg_lib_include, geant4_include]
                                   , dependencies       : neutron_bg_lib_deps
                                   , install            : true
                                   )

neutron_bg_exe = executable( 'neutron_bg'
                               , ['neutron_bg.cc']
                               , include_directories: [neutron_bg_lib_include, nain4_include, geant4_include]
                               , dependencies       : neutron_bg_lib_deps
                               , link_with          : neutron_bg_lib
                               , install            : true
                               )

install_headers(neutron_bg_lib_includes)

pkg = import('pkgconfig')
pkg.generate( neutron_bg_lib
            , description: 'Neutron bg of spallations sources'
            , name       : 'neutron_bg'
            , filebase   : 'neutron_bg'
            )
