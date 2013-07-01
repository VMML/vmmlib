Name:		@CPACK_PACKAGE_NAME@
Version:	@VERSION@
Release:	1%{?dist}
Summary:	@CPACK_PACKAGE_DESCRIPTION_SUMMARY@

Group:		Development/Libraries
License:	@CPACK_PACKAGE_LICENSE@
URL:		http://www.equalizergraphics.com/
Source0:	http://www.equalizergraphics.com/downloads/@CMAKE_PROJECT_NAME@-%{version}.tar.gz
#Patch0:		@CMAKE_PROJECT_NAME@-@VERSION@-build-fix.patch
BuildRequires:	cmake

%description
Vector and Matrix Math Library

%package devel
Summary:	Development files for @CMAKE_PROJECT_NAME@
Group:		Development/Libraries
Requires:	%{name} = %{version}-%{release}

%description devel
Development files for the @CMAKE_PROJECT_NAME@.

%prep
%setup -q
%patch0 -p1 -b .build-fix

%build
%cmake

make %{?_smp_mflags}


%install
make install DESTDIR=%{buildroot}
mv %{buildroot}%{_datadir}/%{name}/doc _tmpdoc/


%post -p /sbin/ldconfig


%postun -p /sbin/ldconfig


%files
%doc _tmpdoc/*
%{_bindir}/*
%{_libdir}/lib*.so.*
%{_datadir}/%{name}

%files devel
%{_includedir}/*
%{_libdir}/lib*.so
%{_libdir}/pkgconfig/*.pc

%changelog
- Initial Release
