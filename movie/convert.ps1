if ($args.Length -ne 1) {
    Write-Host "Usage: $($MyInvocation.MyCommand.Name) <name>"
    exit 1
}

$name = $args[0]
New-Item -ItemType Directory -Force -Path "movie\$name"
Write-Host ""
Write-Host "Converting movie\$name.gif to movie\$name\$name.png..."
magick "movie\$name.gif" "movie\$name\$name.png"