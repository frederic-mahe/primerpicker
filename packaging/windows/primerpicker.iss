; Windows installer of the graphical program (Inno Setup 6): the folder
; made by packaging/windows/deploy-gui.sh, installed for the current
; user (no administrator rights needed; the first page offers an
; installation for all users), with a Start menu entry, an optional
; desktop icon and an uninstaller. Unsigned, like the other packages.
;
; usage: ISCC.exe /DAppVersion=X.Y.Z /DAppFolder=FOLDER /DOutputDir=DIR
;        /DOutputName=NAME packaging\windows\primerpicker.iss
; (FOLDER and DIR: absolute paths)

#ifndef AppVersion
  #error AppVersion is not defined
#endif
#ifndef AppFolder
  #error AppFolder is not defined
#endif
#ifndef OutputDir
  #error OutputDir is not defined
#endif
#ifndef OutputName
  #error OutputName is not defined
#endif

[Setup]
; the identity of the program for Windows (installations, updates and
; uninstallation): never change it
AppId={{0388B3F1-61EA-4100-B406-DF25A31B2BEE}
AppName=PrimerPicker
AppVersion={#AppVersion}
AppVerName=PrimerPicker {#AppVersion}
AppPublisher=Frédéric Mahé
AppPublisherURL=https://github.com/frederic-mahe/primerpicker
AppSupportURL=https://github.com/frederic-mahe/primerpicker/issues
AppUpdatesURL=https://github.com/frederic-mahe/primerpicker/releases
VersionInfoVersion={#AppVersion}
DefaultDirName={autopf}\PrimerPicker
DefaultGroupName=PrimerPicker
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
; 64-bit Windows 10 (1809) or later, as Qt 6
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0.17763
LicenseFile={#AppFolder}\LICENSE
SetupIconFile=..\icons\primerpicker.ico
UninstallDisplayIcon={app}\primerpicker-gui.exe
UninstallDisplayName=PrimerPicker {#AppVersion}
OutputDir={#OutputDir}
OutputBaseFilename={#OutputName}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "{#AppFolder}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\PrimerPicker"; Filename: "{app}\primerpicker-gui.exe"
Name: "{autodesktop}\PrimerPicker"; Filename: "{app}\primerpicker-gui.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\primerpicker-gui.exe"; Description: "{cm:LaunchProgram,PrimerPicker}"; Flags: nowait postinstall skipifsilent
