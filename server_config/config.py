import os
import platform
import shutil

print("Configuration du serveur en cours...")
# Détection du système d'exploitation
is_windows = platform.system() == "Windows"

# Chemin vers le fichier config.toml dans le même répertoire que le script
source_config_path = os.path.join(os.path.dirname(__file__), "config.toml")

# Chemin vers le dossier .streamlit en fonction du système d'exploitation
if is_windows:
    streamlit_folder = os.path.join(os.getenv("USERPROFILE"), ".streamlit")
else:
    streamlit_folder = os.path.expanduser("~/.streamlit")

# Chemin vers le fichier config.toml dans le dossier .streamlit
destination_config_path = os.path.join(streamlit_folder, "config.toml")

# Supprimer le fichier config.toml existant dans le dossier .streamlit s'il existe
if os.path.exists(destination_config_path):
    os.remove(destination_config_path)
    print("Le fichier config.toml existant a été supprimé.")

    # Copier le fichier de configuration
    shutil.copy(source_config_path, destination_config_path)

    print("Le nouveau fichier config.toml a été copié dans le dossier .streamlit.")
else:
    try:
        shutil.copy(source_config_path, destination_config_path)
        print("Configuration effectuée avec succès!")
    except Exception as e:
        print(
            "Un problème a eu lieu lors de la modification du fichier de configuration de streamlit. Il est probable que l'application ne se rafraichisse pas lorsqu'un changement est apporté."
        )
        print(e)
