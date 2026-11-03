plugins {
    id("com.android.application")
}

// CI sets GITHUB_RUN_NUMBER so every build installs as an update over the last one.
val buildNumber = System.getenv("GITHUB_RUN_NUMBER")?.toIntOrNull() ?: 1

android {
    namespace = "com.nicksard.appsleep"
    compileSdk = 35

    defaultConfig {
        applicationId = "com.nicksard.appsleep"
        minSdk = 29
        targetSdk = 35
        versionCode = buildNumber
        versionName = "1.0.$buildNumber"
    }

    // A fixed, non-secret key used only for sideloading. Keeping it in the repo
    // means every CI build has the same signature, so new builds install as
    // updates instead of failing with "App not installed".
    signingConfigs {
        create("sideload") {
            storeFile = file("appsleep-debug.keystore")
            storePassword = "android"
            keyAlias = "androiddebugkey"
            keyPassword = "android"
        }
    }

    buildTypes {
        getByName("release") {
            isMinifyEnabled = false
            signingConfig = signingConfigs.getByName("sideload")
        }
        getByName("debug") {
            signingConfig = signingConfigs.getByName("sideload")
        }
    }

    compileOptions {
        sourceCompatibility = JavaVersion.VERSION_17
        targetCompatibility = JavaVersion.VERSION_17
    }
}
