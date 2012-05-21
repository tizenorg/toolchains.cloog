# 
# Do NOT Edit the Auto-generated Part!
# Generated by: spectacle version 0.22
# 
# >> macros
# << macros

Name:       cloog
Summary:    The Chunky Loop Generator
Version:    0.15.9
Release:    1
Group:      System/Libraries
License:    GPLv2+
URL:        http://www.cloog.org
Source0:    %{name}-ppl-%{version}.tar.gz
BuildRequires:  ppl-devel >= 0.10.2
BuildRequires:  gmp-devel >= 4.1.3
BuildRequires:  libtool


%description
CLooG is a software which generates loops for scanning Z-polyhedra. That is,
CLooG finds the code or pseudo-code where each integral point of one or more
parametrized polyhedron or parametrized polyhedra union is reached. CLooG is
designed to avoid control overhead and to produce a very efficient code.



%package ppl-devel
Summary:    Development tools for the ppl based version of Chunky Loop Generator
Group:      Development/Libraries
Requires:   %{name} = %{version}-%{release}
Requires:   %{name}-ppl = %{version}-%{release}
Requires:   ppl-devel >= 0.10, gmp-devel >= 4.1.3

%description ppl-devel
The header files and dynamic shared libraries of the Chunky Loop Generator.


%package ppl
Summary:    Parma Polyhedra Library backend (ppl) based version of the Cloog binaries
Group:      Development/Libraries
Requires:   %{name} = %{version}-%{release}
Requires(post): /sbin/ldconfig
Requires(postun): /sbin/ldconfig

%description ppl
The dynamic shared libraries of the Chunky Loop Generator


%prep
%setup -q -n %{name}-ppl-%{version}

# >> setup
# << setup

%build
# >> build pre
# << build pre

%configure --disable-static \
    --with-ppl

make %{?jobs:-j%jobs}

# >> build post
# << build post
%install
rm -rf %{buildroot}
# >> install pre
# << install pre
%make_install

# >> install post
rm -rf %{buildroot}/%{_infodir}
# << install post

%post ppl -p /sbin/ldconfig

%postun ppl -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
# >> files
# << files


%files ppl-devel
%defattr(-,root,root,-)
# >> files ppl-devel
%{_includedir}/cloog
%{_libdir}/libcloog.so
# << files ppl-devel

%files ppl
%defattr(-,root,root,-)
# >> files ppl
%{_bindir}/cloog
%{_libdir}/libcloog.so.*
# << files ppl

