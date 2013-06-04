package Version_info;
use strict;
use warnings;
  
##############################################################################  
#
# General information for the triPOD software.
# Change log is located below.
#
##############################################################################    

sub new {
    my $class = shift;
    my $self = {
        _name => 'triPOD',
        _version => 'version 1.2',
        _author => 'Joseph D Baugher',
        _author_email => 'jbaughe2(at)jhmi.edu',
        _url => 'https://github.com/jdbaugher/tripod',
        _download_url => 'git://github.com/jdbaugher/tripod.git',
        _webapp_url => 'http://pevsnerlab.kennedykrieger.org/tripod',
        _reference => 'Baugher JD, Baugher BD, Shirley MD, Pevsner J.,
Sensitive and specific detection of mosaic chromosomal 
abnormalities using the Parent-of-Origin-based Detection
(POD) method. BMC Genomics. 2013 May 31; 14:367.
doi:10.1186/1471-2164-14-367.',
        _ref_url => 'http://www.biomedcentral.com/1471-2164/14/367/'
    };
    bless $self, $class;
    return $self
}

sub get_name {
    my $self = shift;
    return $self->{_name}
}

sub get_version {
    my $self = shift;
    return $self->{_version}
}

sub get_author {
    my $self = shift;
    return $self->{_author}
}

sub get_author_email {
    my $self = shift;
    return $self->{_author_email}
}

sub get_url {
    my $self = shift;
    return $self->{_url}
}

sub get_download_url {
    my $self = shift;
    return $self->{_download_url}
}

sub get_webapp_url {
    my $self = shift;
    return $self->{_webapp_url}
}

sub get_reference {
    my $self = shift;
    return $self->{_reference}
}

sub get_ref_url {
    my $self = shift;
    return $self->{_ref_url}
}


1;

##############################################################################
#
# triPOD Change log 
#
##############################################################################

#---------------
# Version 1.2 - 06/04/2013 - minor changes
#---------------
# Modified the time at which the graphics script is kicked off by the 
#    Perl script. Removed PID tracking and use of two temporary files. If
#    graphics are requested, the Perl script will now finish its analysis,
#    kick off the graphics script, and wait for the graphics to finish.
#    This produces a little longer runtime, but adds more stability across
#    different OSs.
# Added a sample data file to the GitHub repository.

#---------------
# Version 1.1 - 06/01/2013 - minor changes
#--------------- 
# Modified system calls for the PID to work with the darwin core OS.

#---------------
# Version 1.0 - 05/31/2013 - Original publication version
#--------------- 