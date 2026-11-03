# AppSleep

A small Android app (built for a Google Pixel, works on Android 10+) that puts chosen
apps "to sleep" during time windows you set — e.g. Instagram and YouTube from
10 PM to 7 AM on weekdays.

While an app is asleep:

- **Opening it sends you back to the home screen** and shows "*App* is asleep — it wakes
  up at 7:00 AM". If a window starts while you're already in the app, you're kicked out
  within about 30 seconds.
- **Its notifications are cleared** (optional, per schedule).

The app collects nothing and has no internet permission.

## Install on your phone

1. On the Pixel, open the **AppSleep (latest build)** release on this repo's
   GitHub *Releases* page and tap **AppSleep.apk**. (Each push to this folder
   rebuilds it; new builds install as updates.)
2. Open the downloaded file. Allow your browser to "install unknown apps" if asked,
   and tap through the Play Protect prompt ("Install anyway") — it's warning because the
   app didn't come from the Play Store.
3. Open **AppSleep** and follow the card at the top:
   - **Accessibility settings → AppSleep → On.** This is how it notices which app is open.
   - If the switch is greyed out with "Restricted setting": go to
     *Settings → Apps → AppSleep*, tap **⋮** (top right) → **Allow restricted settings**,
     then try again. Android does this for every sideloaded app that asks for accessibility.
   - Optional: **Notification access → AppSleep → On** to clear sleeping apps' notifications.
4. Tap **+ Add schedule**, pick the times, days and apps, and save.

## How it works

Android doesn't let one app freeze or force-stop another (only device-owner/MDM apps can).
So AppSleep does what Play Store app blockers do: an Accessibility service watches for
window changes, and when a sleeping app comes to the front it presses Home and shows a
notice. A notification listener dismisses the app's notifications.

| File | What it does |
|---|---|
| `Schedule.java` | A sleep window (times, days, apps) and the overnight-aware "is it active now?" check |
| `ScheduleStore.java` | Saves schedules as JSON in SharedPreferences |
| `SleepGuardService.java` | Accessibility service that sends you home from sleeping apps |
| `NotificationSilencer.java` | Clears sleeping apps' notifications |
| `MainActivity.java` | Setup status + list of schedules |
| `EditScheduleActivity.java` | Name, times, days, notification option and app picker |
| `AsleepActivity.java` | The "this app is asleep" screen |

## Building it yourself

Open `android-appsleep/` in Android Studio, or run `./gradlew assembleRelease` with the
Android SDK installed. The APK is signed with `app/appsleep-debug.keystore`, a throwaway
key kept in the repo on purpose so every build has the same signature.
