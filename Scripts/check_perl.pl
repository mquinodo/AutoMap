
#!/usr/bin/perl
use ExtUtils::Installed;
my $InstalledModules = ExtUtils::Installed->new(); 
foreach my $module ($InstalledModules->modules()) {
	my $version = $InstalledModules->version($module) || "???";
	print "$module — $versionn" . "\n";
}

